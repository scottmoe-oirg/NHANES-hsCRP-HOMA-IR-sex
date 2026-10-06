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


MARKERS = {
    "ALT": "LBXSATSI",
    "AST": "LBXSASSI",
    "GGT": "LBXSGTSI",
}


def load_config(path: Path):
    with open(path, "rb") as f:
        return tomllib.load(f)


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
        "--hepatic-dir",
        default=str(
            ROOT / "data" / "raw_supplemental"
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

    n_locked = int(domain.sum())

    if n_locked != 3079:
        raise RuntimeError(
            f"Expected locked Study #2 domain n=3079; "
            f"got {n_locked}"
        )

    # --------------------------------------------------------
    # Load hepatic biomarkers
    # --------------------------------------------------------

    hepatic_dir = Path(
        args.hepatic_dir
    )

    files = [
        hepatic_dir / "BIOPRO_I.xpt",
        hepatic_dir / "P_BIOPRO.xpt",
    ]

    frames = []

    needed = [
        "SEQN",
        "LBXSATSI",
        "LBXSASSI",
        "LBXSGTSI",
    ]

    for path in files:

        if not path.exists():
            raise FileNotFoundError(
                f"Missing file: {path}"
            )

        x = pd.read_sas(
            path,
            format="xport",
        )

        missing = [
            col
            for col in needed
            if col not in x.columns
        ]

        if missing:
            raise RuntimeError(
                f"{path.name} missing columns: {missing}"
            )

        frames.append(
            x[needed].copy()
        )

    hepatic = pd.concat(
        frames,
        ignore_index=True,
    )

    if hepatic["SEQN"].duplicated().any():
        raise RuntimeError(
            "Duplicate SEQN in combined hepatic data."
        )

    # --------------------------------------------------------
    # Merge
    # --------------------------------------------------------

    original_rows = len(df)

    df = df.merge(
        hepatic,
        on="SEQN",
        how="left",
        validate="one_to_one",
        sort=False,
    )

    if len(df) != original_rows:
        raise RuntimeError(
            "Hepatic merge changed dataset row count."
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
    # Construct log versions
    # --------------------------------------------------------

    for name, var in MARKERS.items():

        logvar = f"log_{name}"

        analytic[logvar] = np.where(
            analytic[var] > 0,
            np.log(analytic[var]),
            np.nan,
        )

    # --------------------------------------------------------
    # Header
    # --------------------------------------------------------

    print()
    print(
        "STAGE 2E — HEPATIC/METABOLIC BIOCHEMISTRY QC #1"
    )
    print("=" * 105)

    print()
    print(
        f"Locked Study #2 analytic domain: {len(analytic):,}"
    )

    # --------------------------------------------------------
    # Coverage
    # --------------------------------------------------------

    print()
    print("BIOMARKER AVAILABILITY")
    print("-" * 105)

    for name, var in MARKERS.items():

        available = int(
            analytic[var].notna().sum()
        )

        missing = (
            len(analytic) - available
        )

        pct = (
            100
            * available
            / len(analytic)
        )

        print(
            f"{name:>3}: "
            f"available={available:>4,}  "
            f"missing={missing:>3,}  "
            f"coverage={pct:6.2f}%"
        )

    # --------------------------------------------------------
    # Distributions
    # --------------------------------------------------------

    percentiles = [
        .01,
        .025,
        .05,
        .10,
        .25,
        .50,
        .75,
        .90,
        .95,
        .975,
        .99,
    ]

    for name, var in MARKERS.items():

        logvar = f"log_{name}"

        print()
        print(
            f"{name} DISTRIBUTION — RAW"
        )
        print("-" * 105)

        print(
            analytic[var]
            .describe(
                percentiles=percentiles
            )
            .to_string()
        )

        print()
        print(
            f"{name} DISTRIBUTION — LOG"
        )
        print("-" * 105)

        print(
            analytic[logvar]
            .describe(
                percentiles=percentiles
            )
            .to_string()
        )

        print()
        print(
            f"{name}: RAW VS LOG SHAPE"
        )
        print("-" * 105)

        print(
            f"Raw skewness:  "
            f"{analytic[var].skew():.6f}"
        )

        print(
            f"Log skewness:  "
            f"{analytic[logvar].skew():.6f}"
        )

        print(
            f"Raw kurtosis:  "
            f"{analytic[var].kurt():.6f}"
        )

        print(
            f"Log kurtosis:  "
            f"{analytic[logvar].kurt():.6f}"
        )

    # --------------------------------------------------------
    # Availability by sex
    # --------------------------------------------------------

    print()
    print(
        "AVAILABILITY BY SEX — UNWEIGHTED QC"
    )
    print("-" * 105)

    for name, var in MARKERS.items():

        print()
        print(name)

        tab = (
            analytic
            .groupby(
                "sex",
                observed=True,
            )
            .agg(
                n=("SEQN", "size"),
                available=(
                    var,
                    lambda x: x.notna().sum(),
                ),
                missing=(
                    var,
                    lambda x: x.isna().sum(),
                ),
            )
        )

        tab["coverage_pct"] = (
            100
            * tab["available"]
            / tab["n"]
        )

        print(
            tab.to_string()
        )

    # --------------------------------------------------------
    # Availability by period
    # --------------------------------------------------------

    print()
    print(
        "AVAILABILITY BY PERIOD — UNWEIGHTED QC"
    )
    print("-" * 105)

    for name, var in MARKERS.items():

        print()
        print(name)

        tab = (
            analytic
            .groupby(
                "period",
                observed=True,
            )
            .agg(
                n=("SEQN", "size"),
                available=(
                    var,
                    lambda x: x.notna().sum(),
                ),
                missing=(
                    var,
                    lambda x: x.isna().sum(),
                ),
            )
        )

        tab["coverage_pct"] = (
            100
            * tab["available"]
            / tab["n"]
        )

        print(
            tab.to_string()
        )

    # --------------------------------------------------------
    # Distribution by sex
    # --------------------------------------------------------

    print()
    print(
        "RAW BIOMARKER DISTRIBUTIONS BY SEX — "
        "UNWEIGHTED QC"
    )
    print("-" * 105)

    for name, var in MARKERS.items():

        print()
        print(name)

        tab = (
            analytic
            .groupby(
                "sex",
                observed=True,
            )[var]
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
            tab.to_string()
        )

    # --------------------------------------------------------
    # Distribution by period
    # --------------------------------------------------------

    print()
    print(
        "RAW BIOMARKER DISTRIBUTIONS BY PERIOD — "
        "UNWEIGHTED QC"
    )
    print("-" * 105)

    for name, var in MARKERS.items():

        print()
        print(name)

        tab = (
            analytic
            .groupby(
                "period",
                observed=True,
            )[var]
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
            tab.to_string()
        )

    # --------------------------------------------------------
    # Joint availability
    # --------------------------------------------------------

    all_complete = np.ones(
        len(analytic),
        dtype=bool,
    )

    for var in MARKERS.values():
        all_complete &= (
            analytic[var]
            .notna()
            .to_numpy()
        )

    print()
    print(
        "JOINT AVAILABILITY"
    )
    print("-" * 105)

    print(
        "Complete ALT + AST + GGT: "
        f"{int(all_complete.sum()):,}"
    )

    print(
        "Missing at least one:      "
        f"{int((~all_complete).sum()):,}"
    )

    print()
    print("=" * 105)

    print(
        "QC #1 complete. "
        "No hsCRP outcome regression has been fitted."
    )

    print()


if __name__ == "__main__":
    main()
