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


# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

def load_config(path: Path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def egfr_2021_ckd_epi(scr, age, is_male):
    """
    2021 CKD-EPI creatinine equation.

    scr:
        Serum creatinine in mg/dL.

    age:
        Age in years.

    is_male:
        1 = male
        0 = female

    Equation:
        eGFR = 142
               * min(Scr/kappa, 1)^alpha
               * max(Scr/kappa, 1)^(-1.200)
               * 0.9938^Age
               * 1.012 if female

    Sex-specific constants:
        Female: kappa = 0.7, alpha = -0.241
        Male:   kappa = 0.9, alpha = -0.302
    """

    if pd.isna(scr) or pd.isna(age) or pd.isna(is_male):
        return np.nan

    if scr <= 0:
        return np.nan

    if is_male == 1:
        kappa = 0.9
        alpha = -0.302
        sex_factor = 1.0

    elif is_male == 0:
        kappa = 0.7
        alpha = -0.241
        sex_factor = 1.012

    else:
        return np.nan

    ratio = scr / kappa

    return (
        142.0
        * min(ratio, 1.0) ** alpha
        * max(ratio, 1.0) ** (-1.200)
        * (0.9938 ** age)
        * sex_factor
    )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        default=str(
            ROOT
            / "config"
            / "analysis_config.toml"
        ),
    )

    parser.add_argument(
        "--data-dir",
        required=True,
    )

    parser.add_argument(
        "--renal-dir",
        default=str(
            ROOT
            / "data"
            / "raw_supplemental"
        ),
    )

    args = parser.parse_args()

    config = load_config(
        Path(args.config)
    )

    # --------------------------------------------------------
    # Reproduce locked Study #2 dataset
    # --------------------------------------------------------

    df, domain, flow = build_combined_dataset(
        Path(args.data_dir),
        config,
    )

    n_domain = int(domain.sum())

    if n_domain != 3079:
        raise RuntimeError(
            f"Expected locked Study #2 domain n=3079; "
            f"got {n_domain}"
        )

    # --------------------------------------------------------
    # Load serum creatinine
    # --------------------------------------------------------

    renal_dir = Path(args.renal_dir)

    files = [
        (
            "2015-16",
            renal_dir / "BIOPRO_I.xpt",
        ),
        (
            "2017-Mar2020",
            renal_dir / "P_BIOPRO.xpt",
        ),
    ]

    frames = []

    for period, path in files:

        if not path.exists():
            raise FileNotFoundError(
                f"Missing renal file: {path}"
            )

        x = pd.read_sas(
            path,
            format="xport",
        )

        needed = [
            "SEQN",
            "LBXSCR",
            "LBDSCRSI",
        ]

        missing = [
            col
            for col in needed
            if col not in x.columns
        ]

        if missing:
            raise RuntimeError(
                f"{path.name} missing columns: "
                f"{missing}"
            )

        x = x[needed].copy()

        x["renal_period"] = period

        frames.append(x)

    renal = pd.concat(
        frames,
        ignore_index=True,
    )

    if renal["SEQN"].duplicated().any():
        dup = int(
            renal["SEQN"]
            .duplicated()
            .sum()
        )

        raise RuntimeError(
            f"Duplicate SEQN in combined renal data: {dup}"
        )

    # --------------------------------------------------------
    # Merge onto validated Study #2 dataset
    # --------------------------------------------------------

    original_rows = len(df)

    df = df.merge(
        renal,
        on="SEQN",
        how="left",
        validate="one_to_one",
        sort=False,
    )

    if len(df) != original_rows:
        raise RuntimeError(
            "Renal merge changed dataset row count."
        )

    # --------------------------------------------------------
    # Calculate 2021 CKD-EPI eGFR
    # --------------------------------------------------------

    df["eGFR_2021"] = df.apply(
        lambda row: egfr_2021_ckd_epi(
            row["LBXSCR"],
            row["AGE"],
            row["is_male"],
        ),
        axis=1,
    )

    analytic = df.loc[
        domain
    ].copy()

    analytic["sex"] = np.where(
        analytic["is_male"] == 1,
        "Male",
        "Female",
    )

    # --------------------------------------------------------
    # Coverage
    # --------------------------------------------------------

    print()
    print(
        "STAGE 2D — RENAL FUNCTION QC #1"
    )
    print("=" * 100)

    print()
    print("LOCKED STUDY #2 DOMAIN")
    print("-" * 100)

    print(
        f"Analytic n: {len(analytic):,}"
    )

    n_scr = int(
        analytic["LBXSCR"]
        .notna()
        .sum()
    )

    n_missing = (
        len(analytic) - n_scr
    )

    print(
        f"Serum creatinine available: "
        f"{n_scr:,}"
    )

    print(
        f"Serum creatinine missing:   "
        f"{n_missing:,}"
    )

    print(
        f"Coverage: "
        f"{100 * n_scr / len(analytic):.2f}%"
    )

    # --------------------------------------------------------
    # Creatinine distribution
    # --------------------------------------------------------

    print()
    print(
        "SERUM CREATININE DISTRIBUTION (mg/dL)"
    )
    print("-" * 100)

    print(
        analytic["LBXSCR"]
        .describe(
            percentiles=[
                0.01,
                0.025,
                0.05,
                0.25,
                0.50,
                0.75,
                0.95,
                0.975,
                0.99,
            ]
        )
        .to_string()
    )

    # --------------------------------------------------------
    # Log-creatinine distribution
    # --------------------------------------------------------

    analytic["log_LBXSCR"] = np.where(
        analytic["LBXSCR"] > 0,
        np.log(analytic["LBXSCR"]),
        np.nan,
    )

    print()
    print(
        "LOG SERUM CREATININE DISTRIBUTION"
    )
    print("-" * 100)

    print(
        analytic["log_LBXSCR"]
        .describe(
            percentiles=[
                0.01,
                0.025,
                0.05,
                0.10,
                0.25,
                0.50,
                0.75,
                0.90,
                0.95,
                0.975,
                0.99,
            ]
        )
        .to_string()
    )

    print()
    print(
        "RAW VS LOG CREATININE SHAPE"
    )
    print("-" * 100)

    print(
        f"Raw LBXSCR skewness:     "
        f"{analytic['LBXSCR'].skew():.6f}"
    )

    print(
        f"Log LBXSCR skewness:     "
        f"{analytic['log_LBXSCR'].skew():.6f}"
    )

    print(
        f"Raw LBXSCR kurtosis:     "
        f"{analytic['LBXSCR'].kurt():.6f}"
    )

    print(
        f"Log LBXSCR kurtosis:     "
        f"{analytic['log_LBXSCR'].kurt():.6f}"
    )

    # --------------------------------------------------------
    # eGFR distribution
    # --------------------------------------------------------

    print()
    print(
        "2021 CKD-EPI CREATININE eGFR DISTRIBUTION"
    )
    print(
        "(mL/min/1.73 m^2)"
    )
    print("-" * 100)

    print(
        analytic["eGFR_2021"]
        .describe(
            percentiles=[
                0.01,
                0.025,
                0.05,
                0.10,
                0.25,
                0.50,
                0.75,
                0.90,
                0.95,
                0.975,
                0.99,
            ]
        )
        .to_string()
    )

    # --------------------------------------------------------
    # Clinically familiar thresholds
    #
    # Descriptive QC only. These thresholds are NOT being used
    # to define the adjustment variable or exclude participants.
    # --------------------------------------------------------

    egfr = analytic[
        "eGFR_2021"
    ].dropna()

    print()
    print(
        "eGFR THRESHOLDS — UNWEIGHTED QC COUNTS"
    )
    print("-" * 100)

    for threshold in [
        90,
        60,
        45,
        30,
        15,
    ]:

        n = int(
            (egfr < threshold)
            .sum()
        )

        pct = (
            100 * n / len(egfr)
        )

        print(
            f"eGFR < {threshold:>2}: "
            f"{n:>4,} "
            f"({pct:6.2f}%)"
        )

    # --------------------------------------------------------
    # Missingness by sex
    # --------------------------------------------------------

    print()
    print(
        "CREATININE AVAILABILITY BY SEX — "
        "UNWEIGHTED QC"
    )
    print("-" * 100)

    sex_qc = (
        analytic
        .groupby(
            "sex",
            observed=True,
        )
        .agg(
            n=("SEQN", "size"),
            creatinine_available=(
                "LBXSCR",
                lambda x: x.notna().sum(),
            ),
            creatinine_missing=(
                "LBXSCR",
                lambda x: x.isna().sum(),
            ),
        )
    )

    sex_qc["coverage_pct"] = (
        100
        * sex_qc["creatinine_available"]
        / sex_qc["n"]
    )

    print(
        sex_qc.to_string()
    )

    # --------------------------------------------------------
    # Missingness by period
    # --------------------------------------------------------

    print()
    print(
        "CREATININE AVAILABILITY BY PERIOD — "
        "UNWEIGHTED QC"
    )
    print("-" * 100)

    period_qc = (
        analytic
        .groupby(
            "period",
            observed=True,
        )
        .agg(
            n=("SEQN", "size"),
            creatinine_available=(
                "LBXSCR",
                lambda x: x.notna().sum(),
            ),
            creatinine_missing=(
                "LBXSCR",
                lambda x: x.isna().sum(),
            ),
        )
    )

    period_qc["coverage_pct"] = (
        100
        * period_qc["creatinine_available"]
        / period_qc["n"]
    )

    print(
        period_qc.to_string()
    )

    # --------------------------------------------------------
    # eGFR by sex
    # --------------------------------------------------------

    print()
    print(
        "eGFR BY SEX — UNWEIGHTED QC"
    )
    print("-" * 100)

    egfr_sex = (
        analytic
        .groupby(
            "sex",
            observed=True,
        )["eGFR_2021"]
        .agg(
            [
                "count",
                "mean",
                "median",
                "min",
                "max",
            ]
        )
    )

    print(
        egfr_sex.to_string()
    )

    # --------------------------------------------------------
    # eGFR <60 by sex and period
    # --------------------------------------------------------

    analytic["eGFR_lt60"] = np.where(
        analytic["eGFR_2021"].notna(),
        analytic["eGFR_2021"] < 60,
        np.nan,
    )

    print()
    print(
        "eGFR <60 BY SEX — UNWEIGHTED QC"
    )
    print("-" * 100)

    print(
        pd.crosstab(
            analytic["sex"],
            analytic["eGFR_lt60"],
            dropna=False,
        ).to_string()
    )

    print()
    print(
        "eGFR <60 BY PERIOD — UNWEIGHTED QC"
    )
    print("-" * 100)

    print(
        pd.crosstab(
            analytic["period"],
            analytic["eGFR_lt60"],
            dropna=False,
        ).to_string()
    )

    # --------------------------------------------------------
    # Extreme / potentially influential renal values
    # --------------------------------------------------------

    print()
    print(
        "PARTICIPANTS WITH eGFR <30"
    )
    print("-" * 100)

    severe = analytic.loc[
        analytic["eGFR_2021"] < 30,
        [
            "SEQN",
            "sex",
            "AGE",
            "period",
            "LBXSCR",
            "eGFR_2021",
            "HOMA_IR",
            "hs_CRP",
        ],
    ].sort_values(
        "eGFR_2021"
    )

    if len(severe):
        print(
            severe.to_string(
                index=False,
                float_format=lambda x: f"{x:.4f}",
            )
        )
    else:
        print("None.")

    print()
    print("=" * 100)

    print(
        "QC #1 complete. "
        "No outcome regression has been fitted."
    )

    print(
        "Renal specification frozen before outcome modeling: "
        "continuous 2021 CKD-EPI creatinine eGFR as the principal "
        "renal sensitivity, with log serum creatinine as a "
        "complementary directly measured renal specification."
    )

    print()


if __name__ == "__main__":
    main()
