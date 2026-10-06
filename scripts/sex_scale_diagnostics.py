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
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data import build_combined_dataset
from src.survey import survey_linear_regression


def load_config(path: Path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def make_log_homa(df):
    """Construct log HOMA-IR without taking logs of nonpositive values."""
    out = df.copy()
    out["log_HOMA_IR"] = np.nan

    positive = out["HOMA_IR"] > 0
    out.loc[positive, "log_HOMA_IR"] = np.log(
        out.loc[positive, "HOMA_IR"]
    )

    return out


def fit_and_get_residuals(df, domain, outcome, exposure, label):
    """
    Fit the validated survey regression, then reconstruct the Patsy
    design matrix to obtain fitted values and ordinary response residuals.

    The residuals are used diagnostically; inference remains based on
    the survey-design covariance estimator.
    """

    covariates = ["LBXGH", "waist", "AGE", "Non_HDL"]
    cov = " + ".join(covariates)

    formula = (
        f"{outcome} ~ {exposure} * is_male "
        f"+ C(period) * ({cov})"
    )

    result = survey_linear_regression(df, formula, domain)

    # Reconstruct exactly the observations/design matrix Patsy uses.
    y, X = patsy.dmatrices(
        formula,
        df,
        return_type="dataframe",
        NA_action="drop",
    )

    # Restrict to the analytic domain.
    domain_aligned = domain.loc[X.index].astype(bool)

    X_domain = X.loc[domain_aligned].copy()
    y_domain = y.loc[domain_aligned].iloc[:, 0].copy()

    # Ensure coefficient order agrees exactly with the fitted model.
    if list(X_domain.columns) != list(result["names"]):
        raise RuntimeError(
            f"Design-matrix mismatch for model: {label}"
        )

    beta = pd.Series(
        np.asarray(result["beta"]),
        index=result["names"],
    )

    fitted = X_domain @ beta
    residual = y_domain - fitted

    diagnostic = df.loc[X_domain.index].copy()
    diagnostic["observed"] = y_domain
    diagnostic["fitted"] = fitted
    diagnostic["residual"] = residual
    diagnostic["sex"] = np.where(
        diagnostic["is_male"] == 1,
        "Male",
        "Female",
    )
    diagnostic["model"] = label
    diagnostic["outcome_scale"] = outcome
    diagnostic["exposure_scale"] = exposure

    return diagnostic


def residual_summary(d):
    rows = []

    for sex, g in d.groupby("sex"):
        r = g["residual"].dropna()

        rows.append(
            {
                "model": d["model"].iloc[0],
                "sex": sex,
                "n": len(r),
                "mean_residual": r.mean(),
                "sd_residual": r.std(ddof=1),
                "skewness": stats.skew(r, bias=False),
                "excess_kurtosis": stats.kurtosis(
                    r,
                    fisher=True,
                    bias=False,
                ),
                "q01": r.quantile(0.01),
                "q05": r.quantile(0.05),
                "median": r.quantile(0.50),
                "q95": r.quantile(0.95),
                "q99": r.quantile(0.99),
                "max_abs_residual": np.abs(r).max(),
            }
        )

    # Also summarize the entire analytic sample.
    r = d["residual"].dropna()

    rows.append(
        {
            "model": d["model"].iloc[0],
            "sex": "Overall",
            "n": len(r),
            "mean_residual": r.mean(),
            "sd_residual": r.std(ddof=1),
            "skewness": stats.skew(r, bias=False),
            "excess_kurtosis": stats.kurtosis(
                r,
                fisher=True,
                bias=False,
            ),
            "q01": r.quantile(0.01),
            "q05": r.quantile(0.05),
            "median": r.quantile(0.50),
            "q95": r.quantile(0.95),
            "q99": r.quantile(0.99),
            "max_abs_residual": np.abs(r).max(),
        }
    )

    return pd.DataFrame(rows)


def save_residual_vs_fitted(d, path):
    fig, ax = plt.subplots(figsize=(8, 6))

    for sex in ["Female", "Male"]:
        g = d[d["sex"] == sex]
        ax.scatter(
            g["fitted"],
            g["residual"],
            s=14,
            alpha=0.35,
            label=sex,
        )

    ax.axhline(0, linewidth=1)
    ax.set_xlabel("Fitted value")
    ax.set_ylabel("Residual")
    ax.set_title(
        f"Residuals vs fitted\n{d['model'].iloc[0]}"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_residual_histogram(d, path):
    fig, ax = plt.subplots(figsize=(8, 6))

    for sex in ["Female", "Male"]:
        r = d.loc[d["sex"] == sex, "residual"]
        ax.hist(
            r,
            bins=50,
            alpha=0.45,
            density=True,
            label=sex,
        )

    ax.set_xlabel("Residual")
    ax.set_ylabel("Density")
    ax.set_title(
        f"Residual distribution\n{d['model'].iloc[0]}"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_qq(d, path, sex):
    r = d.loc[d["sex"] == sex, "residual"].dropna()

    fig, ax = plt.subplots(figsize=(6, 6))
    stats.probplot(r, dist="norm", plot=ax)

    ax.set_title(
        f"Normal Q-Q: {sex}\n{d['model'].iloc[0]}"
    )

    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_residual_vs_exposure(d, path):
    exposure = d["exposure_scale"].iloc[0]

    fig, ax = plt.subplots(figsize=(8, 6))

    for sex in ["Female", "Male"]:
        g = d[d["sex"] == sex]

        ax.scatter(
            g[exposure],
            g["residual"],
            s=14,
            alpha=0.35,
            label=sex,
        )

    ax.axhline(0, linewidth=1)
    ax.set_xlabel(exposure)
    ax.set_ylabel("Residual")
    ax.set_title(
        f"Residuals vs exposure\n{d['model'].iloc[0]}"
    )
    ax.legend()

    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        default=str(ROOT / "config" / "analysis_config.toml"),
    )

    parser.add_argument(
        "--data-dir",
        required=True,
    )

    parser.add_argument(
        "--out-dir",
        default=str(
            ROOT / "output" / "sex_scale_diagnostics"
        ),
    )

    args = parser.parse_args()

    config = load_config(Path(args.config))

    df, domain, flow = build_combined_dataset(
        Path(args.data_dir),
        config,
    )

    analytic_n = int(domain.sum())

    if analytic_n != 3079:
        raise RuntimeError(
            f"Expected analytic n=3079; got {analytic_n}"
        )

    df = make_log_homa(df)

    if not np.isfinite(
        df.loc[domain, "log_HOMA_IR"]
    ).all():
        raise RuntimeError(
            "Non-finite log HOMA-IR in analytic domain."
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

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    all_summaries = []

    print()
    print("STAGE 1B: SCALE DIAGNOSTICS")
    print("=" * 78)
    print(f"Analytic n: {analytic_n:,}")
    print()

    for i, (label, outcome, exposure) in enumerate(
        models,
        start=1,
    ):

        print(f"[{i}/4] {label}")

        d = fit_and_get_residuals(
            df,
            domain,
            outcome,
            exposure,
            label,
        )

        summary = residual_summary(d)
        all_summaries.append(summary)

        stem = f"model_{i}"

        d.to_csv(
            out / f"{stem}_diagnostic_data.csv",
            index=False,
        )

        save_residual_vs_fitted(
            d,
            out / f"{stem}_residual_vs_fitted.png",
        )

        save_residual_histogram(
            d,
            out / f"{stem}_residual_histogram.png",
        )

        save_residual_vs_exposure(
            d,
            out / f"{stem}_residual_vs_exposure.png",
        )

        save_qq(
            d,
            out / f"{stem}_qq_female.png",
            "Female",
        )

        save_qq(
            d,
            out / f"{stem}_qq_male.png",
            "Male",
        )

    summary = pd.concat(
        all_summaries,
        ignore_index=True,
    )

    summary.to_csv(
        out / "residual_summary.csv",
        index=False,
    )

    print()
    print("RESIDUAL SUMMARY")
    print("-" * 78)

    display_cols = [
        "model",
        "sex",
        "n",
        "sd_residual",
        "skewness",
        "excess_kurtosis",
        "q01",
        "q99",
        "max_abs_residual",
    ]

    print(
        summary[display_cols].to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
        )
    )

    print()
    print(f"Diagnostic plots written to: {out}")


if __name__ == "__main__":
    main()
