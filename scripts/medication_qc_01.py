#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "data" / "raw_supplemental"

FILES = {
    "2015-2016": SUPP / "RXQ_RX_I.xpt",
    "2017-March2020": SUPP / "P_RXQ_RX.xpt",
}


print()
print("MEDICATION QC — RAW FILE STRUCTURE")
print("=" * 90)


for period, path in FILES.items():

    print()
    print(period)
    print("-" * 90)

    if not path.exists():
        raise FileNotFoundError(
            f"Missing file: {path}"
        )

    df = pd.read_sas(
        path,
        format="xport",
    )

    print(f"File: {path.name}")
    print(f"Rows: {len(df):,}")
    print(
        f"Unique SEQN: "
        f"{df['SEQN'].nunique():,}"
    )

    print(
        f"Rows minus unique participants: "
        f"{len(df) - df['SEQN'].nunique():,}"
    )

    print()
    print("COLUMNS")
    print("-" * 90)

    for col in df.columns:
        print(col)

    # --------------------------------------------------------
    # RXDUSE
    # --------------------------------------------------------

    if "RXDUSE" in df.columns:

        print()
        print("RXDUSE COUNTS — ROW LEVEL")
        print("-" * 90)

        print(
            df["RXDUSE"]
            .value_counts(
                dropna=False
            )
            .sort_index()
            .to_string()
        )

        print()
        print("RXDUSE COUNTS — UNIQUE PARTICIPANTS")
        print("-" * 90)

        person_rxduse = (
            df[
                ["SEQN", "RXDUSE"]
            ]
            .drop_duplicates()
        )

        print(
            person_rxduse["RXDUSE"]
            .value_counts(
                dropna=False
            )
            .sort_index()
            .to_string()
        )

    # --------------------------------------------------------
    # Number of rows per participant
    # --------------------------------------------------------

    counts = (
        df.groupby("SEQN")
        .size()
    )

    print()
    print("MEDICATION ROWS PER PARTICIPANT")
    print("-" * 90)

    print(
        counts.describe(
            percentiles=[
                0.50,
                0.75,
                0.90,
                0.95,
                0.99,
            ]
        ).to_string()
    )

    print()
    print(
        "Maximum medication rows for one participant: "
        f"{counts.max():,}"
    )

    # --------------------------------------------------------
    # Drug identifiers
    # --------------------------------------------------------

    for col in [
        "RXDDRUG",
        "RXDDRGID",
        "RXDCOUNT",
    ]:

        if col not in df.columns:
            continue

        print()
        print(f"{col} — FIRST NONMISSING VALUES")
        print("-" * 90)

        vals = (
            df[col]
            .dropna()
            .head(20)
        )

        print(
            vals.to_string(
                index=False
            )
        )

    # --------------------------------------------------------
    # Show sample records
    # --------------------------------------------------------

    desired = [
        "SEQN",
        "RXDUSE",
        "RXDDRUG",
        "RXDDRGID",
        "RXDCOUNT",
    ]

    available = [
        col for col in desired
        if col in df.columns
    ]

    print()
    print("FIRST 30 RECORDS")
    print("-" * 90)

    print(
        df[available]
        .head(30)
        .to_string(
            index=False
        )
    )


print()
print("=" * 90)
print("QC complete.")
print("No medication exposure has been defined.")
print("No regression model has been fitted.")
print()
