#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import tomllib

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data import build_combined_dataset
from src.survey import (
    survey_linear_regression,
    joint_wald_test,
    linear_contrast,
)


# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

def load_config(path: Path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def classify_smoking(row):
    """
    Pre-specified smoking-status definition.

    Never:
        SMQ020 = 2

    Former:
        SMQ020 = 1 and SMQ040 = 3

    Current:
        SMQ020 = 1 and SMQ040 in {1, 2}

    Refused / don't know / otherwise unclassifiable:
        missing
    """

    ever = row["SMQ020"]
    current = row["SMQ040"]

    if ever == 2:
        return "Never"

    if ever == 1 and current == 3:
        return "Former"

    if ever == 1 and current in (1, 2):
        return "Current"

    return np.nan


def fit_model(df, domain, label, add_smoking=False):

    base = (
        "log_hs_CRP ~ HOMA_IR * is_male "
        "+ C(period) * (LBXGH + waist + AGE + Non_HDL)"
    )

    if add_smoking:
        formula = base + " + C(smoking_status)"
    else:
        formula = base

    result = survey_linear_regression(
        df,
        formula,
        domain,
    )

    interaction_term = "HOMA_IR:is_male"

    test = joint_wald_test(
        result,
        [interaction_term],
    )

    female = linear_contrast(
        result,
        {"HOMA_IR": 1.0},
    )

    male = linear_contrast(
        result,
        {
            "HOMA_IR": 1.0,
            interaction_term: 1.0,
        },
    )

    difference = linear_contrast(
        result,
        {interaction_term: 1.0},
    )

    return {
        "model": label,
        "n": result["n"],
        "design_df": result["df"],

        "interaction_F": test["F"],
        "interaction_p": test["p"],

        "female_beta": female["beta"],
        "female_se": female["se"],
        "female_ci_low": female["ci_low"],
        "female_ci_high": female["ci_high"],
        "female_p": female["p"],

        "male_beta": male["beta"],
        "male_se": male["se"],
        "male_ci_low": male["ci_low"],
        "male_ci_high": male["ci_high"],
        "male_p": male["p"],

        "male_minus_female": difference["beta"],
        "difference_se": difference["se"],
        "difference_ci_low": difference["ci_low"],
        "difference_ci_high": difference["ci_high"],
        "difference_p": difference["p"],

        "formula": formula,
    }


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        default=str(
            ROOT / "config" / "analysis_config.toml"
        ),
    )

    parser.add_argument(
        "--data-dir",
        required=True,
    )

    parser.add_argument(
        "--smoking-dir",
        default=str(
            ROOT / "data" / "raw_supplemental"
        ),
    )

    parser.add_argument(
        "--out-dir",
        default=str(
            ROOT / "output" / "sex_smoking_sensitivity"
        ),
    )

    args = parser.parse_args()

    config = load_config(
        Path(args.config)
    )

    # --------------------------------------------------------
    # Build validated Study 2 dataset
    # --------------------------------------------------------

    df, domain, flow = build_combined_dataset(
        Path(args.data_dir),
        config,
    )

    if int(domain.sum()) != 3079:
        raise RuntimeError(
            f"Expected primary analytic domain n=3079; "
            f"got {int(domain.sum())}"
        )

    # --------------------------------------------------------
    # Read smoking files
    # --------------------------------------------------------

    smoking_dir = Path(args.smoking_dir)

    files = [
        smoking_dir / "SMQ_I.xpt",
        smoking_dir / "P_SMQ.xpt",
    ]

    frames = []

    for path in files:

        if not path.exists():
            raise FileNotFoundError(
                f"Missing smoking file: {path}"
            )

        smq = pd.read_sas(
            path,
            format="xport",
        )

        needed = [
            "SEQN",
            "SMQ020",
            "SMQ040",
        ]

        missing_columns = [
            x for x in needed
            if x not in smq.columns
        ]

        if missing_columns:
            raise RuntimeError(
                f"{path.name} missing columns: "
                f"{missing_columns}"
            )

        frames.append(
            smq[needed].copy()
        )

    smoking = pd.concat(
        frames,
        ignore_index=True,
    )

    if smoking["SEQN"].duplicated().any():
        raise RuntimeError(
            "Duplicate SEQN detected in combined smoking data."
        )

    smoking["smoking_status"] = smoking.apply(
        classify_smoking,
        axis=1,
    )

    # Explicit categorical ordering.
    # Never is the reference category.
    smoking["smoking_status"] = pd.Categorical(
        smoking["smoking_status"],
        categories=[
            "Never",
            "Former",
            "Current",
        ],
    )

    # --------------------------------------------------------
    # Merge smoking onto validated dataset
    # --------------------------------------------------------

    original_rows = len(df)

    df = df.merge(
        smoking[
            [
                "SEQN",
                "smoking_status",
            ]
        ],
        on="SEQN",
        how="left",
        validate="one_to_one",
        sort=False,
    )

    if len(df) != original_rows:
        raise RuntimeError(
            "Smoking merge changed dataset row count."
        )

    # --------------------------------------------------------
    # Common complete-case domain
    #
    # Both models are fitted to exactly the same participants.
    # Therefore any change in the interaction estimate is due
    # to smoking adjustment rather than sample composition.
    # --------------------------------------------------------

    smoking_complete = (
        df["smoking_status"].notna()
    )

    common_domain = (
        domain
        & smoking_complete
    )

    n_common = int(
        common_domain.sum()
    )

    if n_common != 3075:
        raise RuntimeError(
            f"Expected smoking-complete domain n=3075; "
            f"got {n_common}"
        )

    # --------------------------------------------------------
    # Fit BOTH models on the same n=3075 sample
    # --------------------------------------------------------

    rows = []

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary model, common sample",
            add_smoking=False,
        )
    )

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary + smoking status",
            add_smoking=True,
        )
    )

    results = pd.DataFrame(rows)

    # --------------------------------------------------------
    # Quantify change in focal sex-interaction estimate
    # --------------------------------------------------------

    reference = results.iloc[0]
    smoking_adjusted = results.iloc[1]

    delta_contrast = (
        smoking_adjusted["male_minus_female"]
        - reference["male_minus_female"]
    )

    abs_ref = abs(
        reference["male_minus_female"]
    )

    abs_adjusted = abs(
        smoking_adjusted["male_minus_female"]
    )

    percent_change_magnitude = (
        100
        * (abs_adjusted - abs_ref)
        / abs_ref
    )

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    out = Path(args.out_dir)

    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    results.to_csv(
        out / "smoking_sensitivity.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Display
    # --------------------------------------------------------

    print()
    print("STAGE 2B: SMOKING ADJUSTMENT SENSITIVITY")
    print("=" * 110)

    print(
        f"Locked primary analytic domain: "
        f"{int(domain.sum()):,}"
    )

    print(
        f"Smoking-complete common domain: "
        f"{n_common:,}"
    )

    print(
        f"Excluded for smoking missingness: "
        f"{int(domain.sum()) - n_common:,}"
    )

    print()
    print(
        "Both models below use exactly the same "
        f"n={n_common:,} participants."
    )

    print()
    print("RESULTS")
    print("-" * 110)

    display_cols = [
        "model",
        "n",
        "design_df",
        "interaction_F",
        "interaction_p",
        "female_beta",
        "male_beta",
        "male_minus_female",
        "difference_ci_low",
        "difference_ci_high",
    ]

    print(
        results[
            display_cols
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )

    print()
    print("CHANGE AFTER SMOKING ADJUSTMENT")
    print("-" * 110)

    print(
        "Male-minus-female contrast before smoking adjustment: "
        f"{reference['male_minus_female']:.6f}"
    )

    print(
        "Male-minus-female contrast after smoking adjustment:  "
        f"{smoking_adjusted['male_minus_female']:.6f}"
    )

    print(
        "Absolute change in contrast:                         "
        f"{delta_contrast:.6f}"
    )

    print(
        "Percent change in absolute contrast magnitude:       "
        f"{percent_change_magnitude:.2f}%"
    )

    print()
    print("SEX-SPECIFIC SLOPES")
    print("-" * 110)

    print(
        "Female:"
    )

    print(
        f"  before smoking adjustment = "
        f"{reference['female_beta']:.6f}"
    )

    print(
        f"  after smoking adjustment  = "
        f"{smoking_adjusted['female_beta']:.6f}"
    )

    print()

    print(
        "Male:"
    )

    print(
        f"  before smoking adjustment = "
        f"{reference['male_beta']:.6f}"
    )

    print(
        f"  after smoking adjustment  = "
        f"{smoking_adjusted['male_beta']:.6f}"
    )

    print()
    print("MODELS")
    print("-" * 110)

    for _, row in results.iterrows():
        print()
        print(row["model"])
        print(row["formula"])

    print()
    print(
        "Smoking definition was pre-specified before "
        "fitting this sensitivity model:"
    )

    print(
        "Never = <100 lifetime cigarettes; "
        "Former = >=100 lifetime cigarettes, not currently smoking; "
        "Current = >=100 lifetime cigarettes and smoking every day/some days."
    )

    print()
    print(
        "Interpretation target: change in the HOMA_IR:is_male "
        "coefficient and its confidence interval."
    )

    print(
        "Do not interpret this adjustment comparison as "
        "mediation or causal evidence."
    )

    print()
    print(
        f"Output: {out / 'smoking_sensitivity.csv'}"
    )

    print()


if __name__ == "__main__":
    main()
