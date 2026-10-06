#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import tomllib

import pandas as pd

# Project root: nhanes_metabolic_followup/
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


def weighted_sex_summary(df, domain):
    """
    Unweighted analytic n and survey-weighted population distribution by sex.
    """
    analytic = df.loc[domain].copy()

    rows = []
    total_weight = analytic["analysis_weight"].sum()

    for label, value in [("Female", 0.0), ("Male", 1.0)]:
        mask = analytic["is_male"].eq(value)
        w = analytic.loc[mask, "analysis_weight"]

        rows.append(
            {
                "sex": label,
                "n_unweighted": int(mask.sum()),
                "weighted_population": float(w.sum()),
                "weighted_percent": float(100.0 * w.sum() / total_weight),
            }
        )

    return pd.DataFrame(rows)


def fit_sex_interaction(df, domain):
    """
    Reconnaissance model:
        log(hs-CRP) ~ HOMA-IR * sex
                      + period-specific nuisance relationships for
                        HbA1c, waist, age, and non-HDL cholesterol

    is_male:
        0 = female
        1 = male

    Therefore:
        HOMA_IR coefficient = female HOMA-IR slope
        HOMA_IR:is_male     = male slope - female slope
    """

    covariates = ["LBXGH", "waist", "AGE", "Non_HDL"]
    cov = " + ".join(covariates)

    formula = (
        f"log_hs_CRP ~ HOMA_IR * is_male "
        f"+ C(period) * ({cov})"
    )

    result = survey_linear_regression(df, formula, domain)

    interaction_term = "HOMA_IR:is_male"
    test = joint_wald_test(result, [interaction_term])

    # Female: is_male = 0
    female = linear_contrast(
        result,
        {"HOMA_IR": 1.0},
    )
    female["sex"] = "Female"

    # Male: is_male = 1
    male = linear_contrast(
        result,
        {
            "HOMA_IR": 1.0,
            "HOMA_IR:is_male": 1.0,
        },
    )
    male["sex"] = "Male"

    slopes = pd.DataFrame([female, male])
    slopes["n_model"] = result["n"]
    slopes["design_df"] = result["df"]

    return formula, result, test, slopes


def main():
    parser = argparse.ArgumentParser(
        description="Initial HOMA-IR x sex reconnaissance analysis."
    )
    parser.add_argument(
        "--config",
        default=str(ROOT / "config" / "analysis_config.toml"),
    )
    parser.add_argument(
        "--data-dir",
        default=str(ROOT / "data" / "raw"),
        help="Directory containing NHANES XPT files",
    )
    parser.add_argument(
        "--out-dir",
        default=str(ROOT / "output" / "sex_reconnaissance"),
    )
    args = parser.parse_args()

    config = load_config(Path(args.config))

    # Reconstruct the exact Study #1 combined dataset and analytic domain.
    df, domain, flow = build_combined_dataset(
        Path(args.data_dir),
        config,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Safety check: this analysis was prospectively defined on the
    # corrected Study #1 analytic population.
    analytic_n = int(domain.sum())
    if analytic_n != 3079:
        raise RuntimeError(
            f"Expected analytic n = 3,079; observed n = {analytic_n}."
        )

    sex_summary = weighted_sex_summary(df, domain)

    formula, result, test, slopes = fit_sex_interaction(
        df,
        domain,
    )

    # Save outputs.
    sex_summary.to_csv(
        out_dir / "sex_descriptive_summary.csv",
        index=False,
    )

    slopes.to_csv(
        out_dir / "sex_specific_HOMA_slopes.csv",
        index=False,
    )

    pd.DataFrame([test]).to_csv(
        out_dir / "HOMA_sex_interaction_test.csv",
        index=False,
    )

    with open(out_dir / "model_formula.txt", "w") as f:
        f.write(formula + "\n")

    # Console report.
    print()
    print("HOMA-IR x SEX RECONNAISSANCE")
    print("=" * 72)
    print(f"Analytic n: {analytic_n:,}")
    print(f"Design df: {result['df']}")
    print()
    print("Analytic population by sex:")
    print(
        sex_summary.to_string(
            index=False,
            formatters={
                "weighted_population": "{:,.0f}".format,
                "weighted_percent": "{:.2f}".format,
            },
        )
    )

    print()
    print("Formal HOMA-IR x sex interaction:")
    print(
        f"F({test['df_num']},{test['df_den']}) = "
        f"{test['F']:.6f}, p = {test['p']:.9f}"
    )

    print()
    print("Adjusted HOMA-IR slopes:")
    print(
        slopes[
            ["sex", "beta", "se", "ci_low", "ci_high", "p"]
        ].to_string(index=False)
    )

    print()
    print("Model:")
    print(formula)

    print()
    print(f"Outputs: {out_dir}")


if __name__ == "__main__":
    main()
