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


def fit_sex_model(df, domain, outcome, exposure, label):
    """
    Fit one cell of the pre-specified 2x2 scale-sensitivity matrix.

    Sex coding:
        is_male = 0 -> Female
        is_male = 1 -> Male

    Therefore:
        exposure coefficient = Female slope
        exposure:is_male     = Male slope - Female slope
    """

    covariates = ["LBXGH", "waist", "AGE", "Non_HDL"]
    cov = " + ".join(covariates)

    formula = (
        f"{outcome} ~ {exposure} * is_male "
        f"+ C(period) * ({cov})"
    )

    result = survey_linear_regression(df, formula, domain)

    interaction_term = f"{exposure}:is_male"
    test = joint_wald_test(result, [interaction_term])

    female = linear_contrast(
        result,
        {exposure: 1.0},
    )

    male = linear_contrast(
        result,
        {
            exposure: 1.0,
            interaction_term: 1.0,
        },
    )

    difference = linear_contrast(
        result,
        {interaction_term: 1.0},
    )

    return {
        "model": label,
        "outcome": outcome,
        "exposure": exposure,
        "formula": formula,
        "n": result["n"],
        "design_df": result["df"],
        "interaction_F": test["F"],
        "interaction_df_num": test["df_num"],
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
    }


def main():
    parser = argparse.ArgumentParser(
        description="Stage 1 sex-interaction scale sensitivity analysis."
    )
    parser.add_argument(
        "--config",
        default=str(ROOT / "config" / "analysis_config.toml"),
    )
    parser.add_argument(
        "--data-dir",
        default=str(ROOT / "data_raw"),
        help="Directory containing NHANES XPT files",
    )
    parser.add_argument(
        "--out-dir",
        default=str(ROOT / "output" / "sex_scale_sensitivity"),
    )
    args = parser.parse_args()

    config = load_config(Path(args.config))

    # Reconstruct the validated Study 1 analytic population.
    df, domain, flow = build_combined_dataset(
        Path(args.data_dir),
        config,
    )

    analytic_n = int(domain.sum())
    if analytic_n != 3079:
        raise RuntimeError(
            f"Expected analytic n = 3,079; observed n = {analytic_n}."
        )

    # Construct log HOMA-IR explicitly.
    # Study 1 required positive HOMA-IR, so this should be defined
    # throughout the analytic domain.
    df = df.copy()
    df["log_HOMA_IR"] = np.log(df["HOMA_IR"])

    if not np.isfinite(df.loc[domain, "log_HOMA_IR"]).all():
        raise RuntimeError(
            "Non-finite log_HOMA_IR values detected in analytic domain."
        )

    models = [
        (
            "log hsCRP / raw HOMA-IR",
            "log_hs_CRP",
            "HOMA_IR",
        ),
        (
            "raw hsCRP / raw HOMA-IR",
            "hs_CRP",
            "HOMA_IR",
        ),
        (
            "log hsCRP / log HOMA-IR",
            "log_hs_CRP",
            "log_HOMA_IR",
        ),
        (
            "raw hsCRP / log HOMA-IR",
            "hs_CRP",
            "log_HOMA_IR",
        ),
    ]

    rows = []

    for label, outcome, exposure in models:
        rows.append(
            fit_sex_model(
                df,
                domain,
                outcome,
                exposure,
                label,
            )
        )

    results = pd.DataFrame(rows)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results.to_csv(
        out_dir / "sex_scale_sensitivity_full.csv",
        index=False,
    )

    # Compact interaction matrix/table.
    interaction_table = results[
        [
            "model",
            "n",
            "design_df",
            "interaction_F",
            "interaction_p",
            "male_minus_female",
            "difference_ci_low",
            "difference_ci_high",
        ]
    ].copy()

    interaction_table.to_csv(
        out_dir / "sex_scale_interaction_summary.csv",
        index=False,
    )

    # Sex-specific slopes.
    slope_rows = []

    for _, r in results.iterrows():
        slope_rows.append(
            {
                "model": r["model"],
                "sex": "Female",
                "beta": r["female_beta"],
                "se": r["female_se"],
                "ci_low": r["female_ci_low"],
                "ci_high": r["female_ci_high"],
                "p": r["female_p"],
            }
        )
        slope_rows.append(
            {
                "model": r["model"],
                "sex": "Male",
                "beta": r["male_beta"],
                "se": r["male_se"],
                "ci_low": r["male_ci_low"],
                "ci_high": r["male_ci_high"],
                "p": r["male_p"],
            }
        )

    slopes = pd.DataFrame(slope_rows)

    slopes.to_csv(
        out_dir / "sex_scale_specific_slopes.csv",
        index=False,
    )

    print()
    print("STAGE 1: SEX-INTERACTION SCALE SENSITIVITY")
    print("=" * 78)
    print(f"Analytic n: {analytic_n:,}")
    print()

    for _, r in results.iterrows():
        print(r["model"])
        print("-" * 78)

        print(
            f"Interaction: "
            f"F({int(r['interaction_df_num'])},"
            f"{int(r['design_df'])}) = "
            f"{r['interaction_F']:.6f}, "
            f"p = {r['interaction_p']:.9f}"
        )

        print(
            "Female slope: "
            f"{r['female_beta']:.6f} "
            f"(SE {r['female_se']:.6f}; "
            f"95% CI {r['female_ci_low']:.6f} to "
            f"{r['female_ci_high']:.6f}; "
            f"p = {r['female_p']:.9f})"
        )

        print(
            "Male slope:   "
            f"{r['male_beta']:.6f} "
            f"(SE {r['male_se']:.6f}; "
            f"95% CI {r['male_ci_low']:.6f} to "
            f"{r['male_ci_high']:.6f}; "
            f"p = {r['male_p']:.9f})"
        )

        print(
            "Male - Female: "
            f"{r['male_minus_female']:.6f} "
            f"(95% CI {r['difference_ci_low']:.6f} to "
            f"{r['difference_ci_high']:.6f})"
        )

        print()

    print("Interaction summary")
    print("-" * 78)

    print(
        interaction_table[
            [
                "model",
                "interaction_F",
                "interaction_p",
            ]
        ].to_string(index=False)
    )

    print()
    print(f"Outputs: {out_dir}")


if __name__ == "__main__":
    main()
