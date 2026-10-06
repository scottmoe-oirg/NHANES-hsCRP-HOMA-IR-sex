#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import tomllib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy

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


def make_log_homa(df):
    out = df.copy()
    out["log_HOMA_IR"] = np.nan

    positive = out["HOMA_IR"] > 0
    out.loc[positive, "log_HOMA_IR"] = np.log(
        out.loc[positive, "HOMA_IR"]
    )

    return out


def fit_model(df, domain, exposure):
    covariates = ["LBXGH", "waist", "AGE", "Non_HDL"]
    cov = " + ".join(covariates)

    formula = (
        f"log_hs_CRP ~ {exposure} * is_male "
        f"+ C(period) * ({cov})"
    )

    result = survey_linear_regression(
        df,
        formula,
        domain,
    )

    return formula, result


def diagnostic_frame(df, domain, formula, result):
    y, X = patsy.dmatrices(
        formula,
        df,
        return_type="dataframe",
        NA_action="drop",
    )

    domain_aligned = domain.loc[X.index].astype(bool)

    Xd = X.loc[domain_aligned].copy()
    yd = y.loc[domain_aligned].iloc[:, 0].copy()

    if list(Xd.columns) != list(result["names"]):
        raise RuntimeError("Design matrix mismatch.")

    beta = pd.Series(
        np.asarray(result["beta"]),
        index=result["names"],
    )

    fitted = Xd @ beta
    residual = yd - fitted

    d = df.loc[Xd.index].copy()

    d["fitted"] = fitted
    d["residual"] = residual

    d["sex"] = np.where(
        d["is_male"] == 1,
        "Male",
        "Female",
    )

    return d


def make_quantile_bins(d, n_bins=10):
    """
    Construct HOMA-IR quantile bins using the pooled analytic sample.

    These are descriptive bins only. They do not redefine the
    regression model.
    """

    d = d.copy()

    d["homa_bin"] = pd.qcut(
        d["HOMA_IR"],
        q=n_bins,
        duplicates="drop",
    )

    return d


def summarize_bins(d):
    rows = []

    for (sex, interval), g in d.groupby(
        ["sex", "homa_bin"],
        observed=True,
    ):
        rows.append(
            {
                "sex": sex,
                "bin": str(interval),
                "n": len(g),
                "HOMA_IR_mean": g["HOMA_IR"].mean(),
                "HOMA_IR_median": g["HOMA_IR"].median(),
                "residual_mean": g["residual"].mean(),
                "residual_sd": g["residual"].std(ddof=1),
                "residual_se": (
                    g["residual"].std(ddof=1)
                    / np.sqrt(len(g))
                ),
            }
        )

    return pd.DataFrame(rows)


def plot_binned_residuals(summary, path, title):
    fig, ax = plt.subplots(figsize=(8, 6))

    for sex in ["Female", "Male"]:
        g = summary[
            summary["sex"] == sex
        ].sort_values("HOMA_IR_mean")

        ax.errorbar(
            g["HOMA_IR_mean"],
            g["residual_mean"],
            yerr=1.96 * g["residual_se"],
            marker="o",
            capsize=3,
            label=sex,
        )

    ax.axhline(0, linewidth=1)

    ax.set_xlabel("Mean HOMA-IR within pooled quantile bin")
    ax.set_ylabel("Mean residual")
    ax.set_title(title)
    ax.legend()

    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def interaction_result(df, domain, exposure, label):
    formula, result = fit_model(
        df,
        domain,
        exposure,
    )

    interaction_term = f"{exposure}:is_male"

    test = joint_wald_test(
        result,
        [interaction_term],
    )

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
        "analysis": label,
        "exposure": exposure,
        "n": result["n"],
        "design_df": result["df"],
        "interaction_F": test["F"],
        "interaction_p": test["p"],
        "female_beta": female["beta"],
        "female_ci_low": female["ci_low"],
        "female_ci_high": female["ci_high"],
        "female_p": female["p"],
        "male_beta": male["beta"],
        "male_ci_low": male["ci_low"],
        "male_ci_high": male["ci_high"],
        "male_p": male["p"],
        "male_minus_female": difference["beta"],
        "difference_ci_low": difference["ci_low"],
        "difference_ci_high": difference["ci_high"],
    }


def upper_tail_thresholds(df, domain):
    x = df.loc[domain, "HOMA_IR"]

    return {
        "Full sample": np.inf,
        "Below 99th percentile": x.quantile(0.99),
        "Below 97.5th percentile": x.quantile(0.975),
        "Below 95th percentile": x.quantile(0.95),
    }


def run_tail_sensitivity(df, domain):
    """
    Diagnostic influence analysis.

    These restrictions are not proposed exclusions from the primary
    analysis. They ask whether the interaction changes when the
    sparsest upper exposure tail is progressively omitted.
    """

    thresholds = upper_tail_thresholds(
        df,
        domain,
    )

    rows = []

    for label, threshold in thresholds.items():

        if np.isinf(threshold):
            restricted_domain = domain.copy()
            cutoff = np.nan
        else:
            restricted_domain = (
                domain
                & (df["HOMA_IR"] <= threshold)
            )
            cutoff = threshold

        for exposure in [
            "HOMA_IR",
            "log_HOMA_IR",
        ]:
            r = interaction_result(
                df,
                restricted_domain,
                exposure,
                label,
            )

            r["HOMA_IR_cutoff"] = cutoff
            r["excluded_n"] = (
                int(domain.sum())
                - int(restricted_domain.sum())
            )

            rows.append(r)

    return pd.DataFrame(rows)


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
            ROOT
            / "output"
            / "sex_exposure_shape_influence"
        ),
    )

    args = parser.parse_args()

    config = load_config(
        Path(args.config)
    )

    df, domain, flow = build_combined_dataset(
        Path(args.data_dir),
        config,
    )

    if int(domain.sum()) != 3079:
        raise RuntimeError(
            f"Expected n=3079; got {int(domain.sum())}"
        )

    df = make_log_homa(df)

    if not np.isfinite(
        df.loc[domain, "log_HOMA_IR"]
    ).all():
        raise RuntimeError(
            "Non-finite log HOMA-IR in analytic domain."
        )

    out = Path(args.out_dir)
    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    #
    # PART A:
    # Binned residual behavior for raw-HOMA model
    #
    formula_raw, result_raw = fit_model(
        df,
        domain,
        "HOMA_IR",
    )

    d_raw = diagnostic_frame(
        df,
        domain,
        formula_raw,
        result_raw,
    )

    d_raw = make_quantile_bins(
        d_raw,
        n_bins=10,
    )

    summary_raw = summarize_bins(
        d_raw,
    )

    summary_raw.to_csv(
        out / "binned_residuals_raw_HOMA.csv",
        index=False,
    )

    plot_binned_residuals(
        summary_raw,
        out / "binned_residuals_raw_HOMA.png",
        "Binned residual means: log hsCRP / raw HOMA-IR",
    )

    #
    # PART B:
    # Binned residual behavior for log-HOMA model
    #
    formula_log, result_log = fit_model(
        df,
        domain,
        "log_HOMA_IR",
    )

    d_log = diagnostic_frame(
        df,
        domain,
        formula_log,
        result_log,
    )

    d_log = make_quantile_bins(
        d_log,
        n_bins=10,
    )

    summary_log = summarize_bins(
        d_log,
    )

    summary_log.to_csv(
        out / "binned_residuals_log_HOMA.csv",
        index=False,
    )

    plot_binned_residuals(
        summary_log,
        out / "binned_residuals_log_HOMA.png",
        "Binned residual means: log hsCRP / log HOMA-IR",
    )

    #
    # PART C:
    # Upper-tail influence sensitivity
    #
    sensitivity = run_tail_sensitivity(
        df,
        domain,
    )

    sensitivity.to_csv(
        out / "upper_tail_sensitivity.csv",
        index=False,
    )

    print()
    print("STAGE 1C: EXPOSURE SHAPE AND UPPER-TAIL INFLUENCE")
    print("=" * 82)
    print(f"Full analytic n: {int(domain.sum()):,}")

    print()
    print("HOMA-IR DISTRIBUTION")
    print("-" * 82)

    x = df.loc[domain, "HOMA_IR"]

    for q in [
        0.50,
        0.75,
        0.90,
        0.95,
        0.975,
        0.99,
        1.00,
    ]:
        print(
            f"{100*q:5.1f}th percentile: "
            f"{x.quantile(q):.6f}"
        )

    print()
    print("UPPER-TAIL SENSITIVITY")
    print("-" * 82)

    cols = [
        "analysis",
        "exposure",
        "n",
        "excluded_n",
        "HOMA_IR_cutoff",
        "interaction_F",
        "interaction_p",
        "female_beta",
        "male_beta",
        "male_minus_female",
    ]

    print(
        sensitivity[cols].to_string(
            index=False,
            float_format=lambda z: f"{z:.6f}",
        )
    )

    print()
    print("BINNED RESIDUAL MEANS: RAW HOMA-IR MODEL")
    print("-" * 82)

    print(
        summary_raw[
            [
                "sex",
                "n",
                "HOMA_IR_mean",
                "residual_mean",
                "residual_se",
            ]
        ].to_string(
            index=False,
            float_format=lambda z: f"{z:.5f}",
        )
    )

    print()
    print("BINNED RESIDUAL MEANS: LOG HOMA-IR MODEL")
    print("-" * 82)

    print(
        summary_log[
            [
                "sex",
                "n",
                "HOMA_IR_mean",
                "residual_mean",
                "residual_se",
            ]
        ].to_string(
            index=False,
            float_format=lambda z: f"{z:.5f}",
        )
    )

    print()
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
