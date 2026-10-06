#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import tomllib

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


def fit_model(df, domain, label):

    formula = (
        "log_hs_CRP ~ HOMA_IR * is_male "
        "+ C(period) * (LBXGH + waist + AGE + Non_HDL)"
    )

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
        "--out-dir",
        default=str(
            ROOT / "output" / "sex_hscrp10_sensitivity"
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

    n_primary = int(domain.sum())

    if n_primary != 3079:
        raise RuntimeError(
            f"Expected primary analytic domain n=3079; "
            f"got {n_primary}"
        )

    # --------------------------------------------------------
    # Frozen hsCRP restriction
    #
    # Prespecified sensitivity:
    #     raw hsCRP <= 10 mg/L
    #
    # The primary outcome remains log(hsCRP).
    # No alternate hsCRP threshold or outcome scale is examined.
    # --------------------------------------------------------

    if df.loc[domain, "hs_CRP"].isna().any():
        raise RuntimeError(
            "Missing hsCRP detected in locked primary domain."
        )

    if (df.loc[domain, "hs_CRP"] <= 0).any():
        raise RuntimeError(
            "Nonpositive hsCRP detected in locked primary domain."
        )

    restricted_domain = (
        domain
        & (df["hs_CRP"] <= 10.0)
    )

    n_restricted = int(
        restricted_domain.sum()
    )

    n_excluded = (
        n_primary - n_restricted
    )

    if n_restricted != 2895:
        raise RuntimeError(
            f"Expected hsCRP<=10 domain n=2895; "
            f"got {n_restricted}"
        )

    if n_excluded != 184:
        raise RuntimeError(
            f"Expected 184 observations with hsCRP>10; "
            f"got {n_excluded}"
        )

    # --------------------------------------------------------
    # Fit primary model on:
    #
    # 1. Full locked domain
    # 2. Frozen hsCRP <=10 mg/L domain
    #
    # Survey design is preserved because observations outside
    # each analytic domain remain in the dataframe and their
    # estimating-equation contributions are score-zeroed by
    # the validated survey regression machinery.
    # --------------------------------------------------------

    rows = []

    rows.append(
        fit_model(
            df,
            domain,
            "Primary model, full domain",
        )
    )

    rows.append(
        fit_model(
            df,
            restricted_domain,
            "Primary model, hsCRP <=10 mg/L",
        )
    )

    results = pd.DataFrame(rows)

    reference = results.iloc[0]
    restricted = results.iloc[1]

    # --------------------------------------------------------
    # Quantify change in focal sex-interaction estimate
    # --------------------------------------------------------

    delta_contrast = (
        restricted["male_minus_female"]
        - reference["male_minus_female"]
    )

    abs_ref = abs(
        reference["male_minus_female"]
    )

    abs_restricted = abs(
        restricted["male_minus_female"]
    )

    percent_change_magnitude = (
        100.0
        * (abs_restricted - abs_ref)
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
        out / "hscrp10_sensitivity.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Display
    # --------------------------------------------------------

    print()
    print("STAGE 2H: hsCRP <=10 mg/L SENSITIVITY")
    print("=" * 110)

    print(
        f"Locked primary analytic domain: "
        f"{n_primary:,}"
    )

    print(
        f"Participants with hsCRP >10 mg/L: "
        f"{n_excluded:,}"
    )

    print(
        f"Restricted hsCRP <=10 mg/L domain: "
        f"{n_restricted:,}"
    )

    print()
    print(
        "Outcome remains log(hsCRP); only the analytic "
        "domain is restricted."
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
    print("CHANGE AFTER hsCRP RESTRICTION")
    print("-" * 110)

    print(
        "Male-minus-female contrast, full domain:       "
        f"{reference['male_minus_female']:.6f}"
    )

    print(
        "Male-minus-female contrast, hsCRP <=10 mg/L:  "
        f"{restricted['male_minus_female']:.6f}"
    )

    print(
        "Restricted - full contrast:                   "
        f"{delta_contrast:.6f}"
    )

    print(
        "Percent change in absolute contrast magnitude: "
        f"{percent_change_magnitude:.2f}%"
    )

    print()
    print("SEX-SPECIFIC SLOPES")
    print("-" * 110)

    print("Female:")

    print(
        f"  full domain       = "
        f"{reference['female_beta']:.6f}"
    )

    print(
        f"  hsCRP <=10 mg/L   = "
        f"{restricted['female_beta']:.6f}"
    )

    print()

    print("Male:")

    print(
        f"  full domain       = "
        f"{reference['male_beta']:.6f}"
    )

    print(
        f"  hsCRP <=10 mg/L   = "
        f"{restricted['male_beta']:.6f}"
    )

    print()
    print("MODELS")
    print("-" * 110)

    for _, row in results.iterrows():

        print()
        print(row["model"])
        print(row["formula"])

    print()
    print("INTERPRETIVE NOTE")
    print("-" * 110)

    print(
        "The hsCRP <=10 mg/L restriction was specified "
        "before examining the Study #2 sex-interaction result "
        "under this restriction."
    )

    print(
        "hsCRP >10 mg/L is not interpreted as synonymous "
        "with acute infection; this is a robustness restriction "
        "for high inflammatory states."
    )

    print()
    print(
        f"Saved: {out / 'hscrp10_sensitivity.csv'}"
    )


if __name__ == "__main__":
    main()
