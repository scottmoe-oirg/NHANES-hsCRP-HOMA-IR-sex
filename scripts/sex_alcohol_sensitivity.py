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


def normalize_xpt_zero(x):
    """
    SAS/XPORT numeric zero may appear as an extremely small
    positive floating-point value (~5.397605e-79).
    """
    if pd.notna(x) and abs(float(x)) < 1e-50:
        return 0.0
    return x


def classify_old(row):
    """
    Frozen 2015-2016 alcohol-status definition.

    Zero:
        ALQ120Q = 0
        Explicitly no drinking occasions in past 12 months.

    Very-low-lifetime:
        ALQ101 = 2 and ALQ110 = 2
        Fewer than 12 drinks in any year AND fewer than
        12 drinks in entire life. Exact past-year volume
        cannot be reconstructed.

    Positive:
        Valid ALQ120Q > 0
        Positive drinking frequency in past 12 months.

    Otherwise:
        missing / genuinely unclassifiable
    """

    q = row["ALQ120Q"]

    if pd.notna(q) and q == 0:
        return "Zero"

    if (
        row["ALQ101"] == 2
        and row["ALQ110"] == 2
    ):
        return "Very-low lifetime"

    if (
        pd.notna(q)
        and q > 0
        and q not in (777, 999)
    ):
        return "Positive"

    return np.nan


def classify_new(row):
    """
    Frozen 2017-March 2020 alcohol-status definition.

    Zero:
        ALQ111 = 2
        OR ALQ121 = 0.

    Positive:
        Valid ALQ121 > 0.

    Otherwise:
        missing / genuinely unclassifiable.
    """

    if row["ALQ111"] == 2:
        return "Zero"

    q = row["ALQ121"]

    if pd.notna(q) and q == 0:
        return "Zero"

    if (
        pd.notna(q)
        and q > 0
        and q not in (77, 99)
    ):
        return "Positive"

    return np.nan


def fit_model(df, domain, label, add_alcohol=False):

    base = (
        "log_hs_CRP ~ HOMA_IR * is_male "
        "+ C(period) * (LBXGH + waist + AGE + Non_HDL)"
    )

    if add_alcohol:
        formula = base + " + C(alcohol_status)"
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
        "--alcohol-dir",
        default=str(
            ROOT / "data" / "raw_supplemental"
        ),
    )

    parser.add_argument(
        "--out-dir",
        default=str(
            ROOT / "output" / "sex_alcohol_sensitivity"
        ),
    )

    args = parser.parse_args()

    config = load_config(
        Path(args.config)
    )

    # --------------------------------------------------------
    # Build validated Study #2 dataset
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
    # Read alcohol files
    # --------------------------------------------------------

    alcohol_dir = Path(args.alcohol_dir)

    old_path = alcohol_dir / "ALQ_I.xpt"
    new_path = alcohol_dir / "P_ALQ.xpt"

    if not old_path.exists():
        raise FileNotFoundError(
            f"Missing alcohol file: {old_path}"
        )

    if not new_path.exists():
        raise FileNotFoundError(
            f"Missing alcohol file: {new_path}"
        )

    old = pd.read_sas(
        old_path,
        format="xport",
    )

    new = pd.read_sas(
        new_path,
        format="xport",
    )

    old_needed = [
        "SEQN",
        "ALQ101",
        "ALQ110",
        "ALQ120Q",
    ]

    new_needed = [
        "SEQN",
        "ALQ111",
        "ALQ121",
    ]

    old_missing = [
        x for x in old_needed
        if x not in old.columns
    ]

    new_missing = [
        x for x in new_needed
        if x not in new.columns
    ]

    if old_missing:
        raise RuntimeError(
            f"{old_path.name} missing columns: "
            f"{old_missing}"
        )

    if new_missing:
        raise RuntimeError(
            f"{new_path.name} missing columns: "
            f"{new_missing}"
        )

    old = old[old_needed].copy()
    new = new[new_needed].copy()

    # --------------------------------------------------------
    # Normalize SAS/XPORT pseudo-zero
    # --------------------------------------------------------

    old["ALQ120Q"] = (
        old["ALQ120Q"]
        .apply(normalize_xpt_zero)
    )

    new["ALQ121"] = (
        new["ALQ121"]
        .apply(normalize_xpt_zero)
    )

    # --------------------------------------------------------
    # Apply FROZEN alcohol definition
    # --------------------------------------------------------

    old["alcohol_status"] = old.apply(
        classify_old,
        axis=1,
    )

    new["alcohol_status"] = new.apply(
        classify_new,
        axis=1,
    )

    alcohol = pd.concat(
        [
            old[
                [
                    "SEQN",
                    "alcohol_status",
                ]
            ],
            new[
                [
                    "SEQN",
                    "alcohol_status",
                ]
            ],
        ],
        ignore_index=True,
    )

    if alcohol["SEQN"].duplicated().any():
        raise RuntimeError(
            "Duplicate SEQN detected in combined alcohol data."
        )

    # --------------------------------------------------------
    # Explicit categorical ordering
    #
    # Zero is reference category.
    # --------------------------------------------------------

    alcohol["alcohol_status"] = pd.Categorical(
        alcohol["alcohol_status"],
        categories=[
            "Zero",
            "Very-low lifetime",
            "Positive",
        ],
    )

    # --------------------------------------------------------
    # Merge alcohol onto validated dataset
    # --------------------------------------------------------

    original_rows = len(df)

    df = df.merge(
        alcohol[
            [
                "SEQN",
                "alcohol_status",
            ]
        ],
        on="SEQN",
        how="left",
        validate="one_to_one",
        sort=False,
    )

    if len(df) != original_rows:
        raise RuntimeError(
            "Alcohol merge changed dataset row count."
        )

    # --------------------------------------------------------
    # Verify frozen QC #4 counts in locked domain
    # --------------------------------------------------------

    locked_status = (
        df.loc[
            domain,
            "alcohol_status",
        ]
        .value_counts(
            dropna=False
        )
    )

    observed_zero = int(
        locked_status.get(
            "Zero",
            0,
        )
    )

    observed_low = int(
        locked_status.get(
            "Very-low lifetime",
            0,
        )
    )

    observed_positive = int(
        locked_status.get(
            "Positive",
            0,
        )
    )

    observed_unknown = (
        int(domain.sum())
        - observed_zero
        - observed_low
        - observed_positive
    )

    expected = {
        "Zero": 555,
        "Very-low lifetime": 163,
        "Positive": 2201,
        "Unknown": 160,
    }

    observed = {
        "Zero": observed_zero,
        "Very-low lifetime": observed_low,
        "Positive": observed_positive,
        "Unknown": observed_unknown,
    }

    if observed != expected:
        raise RuntimeError(
            "Frozen alcohol-state counts do not match QC #4.\n"
            f"Expected: {expected}\n"
            f"Observed: {observed}"
        )

    # --------------------------------------------------------
    # Common complete-case domain
    #
    # Genuine unknown/unclassifiable alcohol status is excluded.
    #
    # Both models are fitted to exactly the same participants.
    # Therefore any change in the interaction estimate is due
    # to alcohol adjustment rather than sample composition.
    # --------------------------------------------------------

    alcohol_complete = (
        df["alcohol_status"].notna()
    )

    common_domain = (
        domain
        & alcohol_complete
    )

    n_common = int(
        common_domain.sum()
    )

    if n_common != 2919:
        raise RuntimeError(
            f"Expected alcohol-classifiable domain n=2919; "
            f"got {n_common}"
        )

    # --------------------------------------------------------
    # Fit BOTH models on same n=2919 sample
    # --------------------------------------------------------

    rows = []

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary model, common sample",
            add_alcohol=False,
        )
    )

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary + alcohol status",
            add_alcohol=True,
        )
    )

    results = pd.DataFrame(rows)

    # --------------------------------------------------------
    # Quantify change in focal sex-interaction estimate
    # --------------------------------------------------------

    reference = results.iloc[0]
    alcohol_adjusted = results.iloc[1]

    delta_contrast = (
        alcohol_adjusted["male_minus_female"]
        - reference["male_minus_female"]
    )

    abs_ref = abs(
        reference["male_minus_female"]
    )

    abs_adjusted = abs(
        alcohol_adjusted["male_minus_female"]
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
        out / "alcohol_sensitivity.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Display
    # --------------------------------------------------------

    print()
    print("STAGE 2G: ALCOHOL ADJUSTMENT SENSITIVITY")
    print("=" * 110)

    print(
        f"Locked primary analytic domain: "
        f"{int(domain.sum()):,}"
    )

    print()
    print("FROZEN ALCOHOL-STATUS COUNTS")
    print("-" * 110)

    for key in [
        "Zero",
        "Very-low lifetime",
        "Positive",
        "Unknown",
    ]:
        print(
            f"{key:20s}: "
            f"{observed[key]:,}"
        )

    print()
    print(
        f"Alcohol-classifiable common domain: "
        f"{n_common:,}"
    )

    print(
        f"Excluded as genuinely unknown/unclassifiable: "
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
    print("CHANGE AFTER ALCOHOL ADJUSTMENT")
    print("-" * 110)

    print(
        "Male-minus-female contrast before alcohol adjustment: "
        f"{reference['male_minus_female']:.6f}"
    )

    print(
        "Male-minus-female contrast after alcohol adjustment:  "
        f"{alcohol_adjusted['male_minus_female']:.6f}"
    )

    print(
        "Absolute change in contrast:                          "
        f"{delta_contrast:.6f}"
    )

    print(
        "Percent change in absolute contrast magnitude:        "
        f"{percent_change_magnitude:.2f}%"
    )

    print()
    print("SEX-SPECIFIC SLOPES")
    print("-" * 110)

    print("Female:")

    print(
        f"  before alcohol adjustment = "
        f"{reference['female_beta']:.6f}"
    )

    print(
        f"  after alcohol adjustment  = "
        f"{alcohol_adjusted['female_beta']:.6f}"
    )

    print()

    print("Male:")

    print(
        f"  before alcohol adjustment = "
        f"{reference['male_beta']:.6f}"
    )

    print(
        f"  after alcohol adjustment  = "
        f"{alcohol_adjusted['male_beta']:.6f}"
    )

    print()
    print("MODELS")
    print("-" * 110)

    for _, row in results.iterrows():
        print()
        print(row["model"])
        print(row["formula"])

    print()
    print("FROZEN ALCOHOL DEFINITION")
    print("-" * 110)

    print(
        "Zero = documented zero past-year alcohol exposure."
    )

    print(
        "Very-low lifetime = 2015-2016 participants reporting "
        "fewer than 12 drinks in any year and fewer than "
        "12 drinks in their lifetime, but without exact "
        "past-year volume."
    )

    print(
        "Positive = documented positive past-year drinking frequency."
    )

    print(
        "Genuinely unknown/unclassifiable alcohol exposure was excluded."
    )

    print()
    print(
        "This categorical specification was frozen before "
        "examining HOMA-IR/hsCRP results for the alcohol sensitivity."
    )

    print()
    print(
        f"Saved: {out / 'alcohol_sensitivity.csv'}"
    )


if __name__ == "__main__":
    main()
