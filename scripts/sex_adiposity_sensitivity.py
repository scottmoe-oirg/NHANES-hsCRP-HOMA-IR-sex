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


def load_config(path: Path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def fit_one(df, domain, label, adiposity_term):
    """
    Fit the locked Study 2 interaction model while varying only
    the adiposity adjustment.

    Locked scientific structure:
        log_hs_CRP ~ HOMA_IR * is_male

    Locked background adjustment:
        LBXGH + AGE + Non_HDL

    Variable component:
        none, BMI, waist, or waist_to_height
    """

    covariates = ["LBXGH", "AGE", "Non_HDL"]

    if adiposity_term is not None:
        covariates.append(adiposity_term)

    cov = " + ".join(covariates)

    formula = (
        f"log_hs_CRP ~ HOMA_IR * is_male "
        f"+ C(period) * ({cov})"
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
        "adiposity_term": (
            "none" if adiposity_term is None else adiposity_term
        ),
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


def print_results(title, results):
    print()
    print(title)
    print("=" * 108)

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
        results[display_cols].to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )


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
            ROOT / "output" / "sex_adiposity_sensitivity"
        ),
    )

    args = parser.parse_args()

    config = load_config(Path(args.config))

    df, domain, flow = build_combined_dataset(
        Path(args.data_dir),
        config,
    )

    if int(domain.sum()) != 3079:
        raise RuntimeError(
            f"Expected primary analytic domain n=3079; "
            f"got {int(domain.sum())}"
        )

    models = [
        ("No adiposity adjustment", None),
        ("BMI adjustment", "BMI"),
        ("Waist adjustment [PRIMARY]", "waist"),
        ("Waist-to-height adjustment", "waist_to_height"),
    ]

    #
    # A. Available-case models
    #
    # survey_linear_regression handles formula-specific missingness.
    #
    available_rows = []

    for label, adiposity_term in models:
        available_rows.append(
            fit_one(
                df,
                domain,
                label,
                adiposity_term,
            )
        )

    available = pd.DataFrame(available_rows)

    #
    # B. Common complete-case domain
    #
    # Require every variable needed by every adiposity model.
    #
    common_vars = [
        "log_hs_CRP",
        "HOMA_IR",
        "is_male",
        "period",
        "LBXGH",
        "AGE",
        "Non_HDL",
        "BMI",
        "waist",
        "waist_to_height",
    ]

    complete = df[common_vars].notna().all(axis=1)

    finite = np.isfinite(
        df[common_vars].select_dtypes(
            include=[np.number]
        )
    ).all(axis=1)

    common_domain = (
        domain
        & complete
        & finite
    )

    common_rows = []

    for label, adiposity_term in models:
        common_rows.append(
            fit_one(
                df,
                common_domain,
                label,
                adiposity_term,
            )
        )

    common = pd.DataFrame(common_rows)

    #
    # Save outputs
    #
    out = Path(args.out_dir)
    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    available.to_csv(
        out / "adiposity_available_case.csv",
        index=False,
    )

    common.to_csv(
        out / "adiposity_common_complete_case.csv",
        index=False,
    )

    #
    # Print useful missingness information
    #
    print()
    print("STAGE 2A: ADIPOSITY REPRESENTATION")
    print("=" * 108)
    print(f"Locked primary analytic domain: {int(domain.sum()):,}")
    print(f"Common complete-case domain:   {int(common_domain.sum()):,}")
    print()

    print("ADIPOSITY VARIABLE AVAILABILITY WITHIN LOCKED DOMAIN")
    print("-" * 108)

    for var in ["BMI", "waist", "height", "waist_to_height"]:
        n_present = int(
            (
                domain
                & df[var].notna()
                & np.isfinite(df[var])
            ).sum()
        )

        print(
            f"{var:18s}: "
            f"{n_present:,} / {int(domain.sum()):,}"
        )

    print_results(
        "A. AVAILABLE-CASE MODELS",
        available,
    )

    print_results(
        "B. COMMON COMPLETE-CASE MODELS",
        common,
    )

    print()
    print("Locked primary specification:")
    print(
        "log_hs_CRP ~ HOMA_IR * is_male "
        "+ C(period) * (LBXGH + waist + AGE + Non_HDL)"
    )

    print()
    print(
        "Interpretation target: change in the HOMA_IR:is_male "
        "coefficient and its CI across adiposity specifications."
    )

    print()
    print(
        "These are adjustment-sensitivity analyses, not formal "
        "mediation analyses."
    )

    print()
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
