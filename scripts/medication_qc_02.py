#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "data" / "raw_supplemental" / "RXQ_DRUG.xpt"


print()
print("MEDICATION QC #2 — DRUG INFORMATION / THERAPEUTIC CLASSIFICATION")
print("=" * 100)


if not PATH.exists():
    raise FileNotFoundError(
        f"Missing file: {PATH}"
    )


df = pd.read_sas(
    PATH,
    format="xport",
)


print()
print("FILE STRUCTURE")
print("-" * 100)

print(f"File: {PATH.name}")
print(f"Rows: {len(df):,}")
print(f"Columns: {len(df.columns):,}")


print()
print("COLUMNS")
print("-" * 100)

for col in df.columns:
    print(col)


# ------------------------------------------------------------------
# General uniqueness / identifier inspection
# ------------------------------------------------------------------

for candidate in [
    "RXDDRGID",
    "RXDDRUG",
]:

    if candidate in df.columns:

        print()
        print(f"{candidate} SUMMARY")
        print("-" * 100)

        print(
            f"Nonmissing: "
            f"{df[candidate].notna().sum():,}"
        )

        print(
            f"Unique: "
            f"{df[candidate].nunique(dropna=True):,}"
        )

        print()
        print("First 25 nonmissing values:")

        print(
            df[candidate]
            .dropna()
            .head(25)
            .to_string(index=False)
        )


# ------------------------------------------------------------------
# Identify columns that look like therapeutic/class variables
# ------------------------------------------------------------------

keywords = [
    "TC",
    "THER",
    "CLASS",
    "CAT",
    "ING",
]

possible_class_cols = []

for col in df.columns:

    upper = col.upper()

    if any(
        key in upper
        for key in keywords
    ):
        possible_class_cols.append(col)


print()
print("POSSIBLE THERAPEUTIC / CLASSIFICATION COLUMNS")
print("-" * 100)

if possible_class_cols:

    for col in possible_class_cols:
        print(col)

else:
    print("None identified by column-name search.")


# ------------------------------------------------------------------
# Show first rows across all columns
# ------------------------------------------------------------------

print()
print("FIRST 20 COMPLETE RECORDS")
print("-" * 100)

with pd.option_context(
    "display.max_columns",
    None,
    "display.width",
    250,
):
    print(
        df.head(20)
        .to_string(index=False)
    )


# ------------------------------------------------------------------
# Search drug-name fields for several known lipid-lowering examples.
#
# This is exploratory QC ONLY.
# It does NOT define the eventual medication exposure.
# ------------------------------------------------------------------

drug_name_cols = [
    col
    for col in df.columns
    if "DRUG" in col.upper()
    or "NAME" in col.upper()
]


examples = [
    "SIMVASTATIN",
    "ATORVASTATIN",
    "ROSUVASTATIN",
    "PRAVASTATIN",
    "GEMFIBROZIL",
    "FENOFIBRATE",
    "EZETIMIBE",
]


print()
print("KNOWN LIPID-LOWERING DRUG EXAMPLES")
print("-" * 100)

if not drug_name_cols:

    print(
        "No obvious drug-name columns identified."
    )

else:

    print(
        "Searching columns:",
        drug_name_cols,
    )

    for drug in examples:

        mask = pd.Series(
            False,
            index=df.index,
        )

        for col in drug_name_cols:

            # Convert bytes/strings safely for QC search.
            values = (
                df[col]
                .astype(str)
                .str.upper()
            )

            mask = (
                mask
                | values.str.contains(
                    drug,
                    na=False,
                    regex=False,
                )
            )

        hits = df.loc[mask]

        print()
        print(
            f"{drug}: "
            f"{len(hits):,} matching row(s)"
        )

        if len(hits):

            with pd.option_context(
                "display.max_columns",
                None,
                "display.width",
                250,
            ):
                print(
                    hits.head(10)
                    .to_string(index=False)
                )


# ------------------------------------------------------------------
# Check whether drug IDs from participant files match this file.
# ------------------------------------------------------------------

print()
print("LINKAGE CHECK WITH PARTICIPANT PRESCRIPTION FILES")
print("-" * 100)

rx_files = [
    ROOT / "data" / "raw_supplemental" / "RXQ_RX_I.xpt",
    ROOT / "data" / "raw_supplemental" / "P_RXQ_RX.xpt",
]

if "RXDDRGID" not in df.columns:

    print(
        "RXDDRGID is not present in RXQ_DRUG; "
        "cannot perform linkage check."
    )

else:

    drug_ids = set(
        df["RXDDRGID"]
        .dropna()
    )

    for path in rx_files:

        rx = pd.read_sas(
            path,
            format="xport",
        )

        ids = (
            rx.loc[
                rx["RXDDRGID"].notna(),
                "RXDDRGID",
            ]
        )

        # Blank byte/string values are technically nonmissing,
        # so exclude them from the meaningful-ID denominator.
        meaningful = ids[
            ids.astype(str).isin(
                [
                    "b''",
                    "''",
                ]
            )
            == False
        ]

        matched = meaningful.isin(
            drug_ids
        )

        print()
        print(path.name)

        print(
            f"Meaningful RXDDRGID records: "
            f"{len(meaningful):,}"
        )

        print(
            f"Matched to RXQ_DRUG: "
            f"{matched.sum():,}"
        )

        if len(meaningful):

            print(
                f"Match percentage: "
                f"{100 * matched.mean():.2f}%"
            )

        unmatched = (
            meaningful[
                ~matched
            ]
            .drop_duplicates()
        )

        print(
            f"Unique unmatched IDs: "
            f"{len(unmatched):,}"
        )

        if len(unmatched):

            print(
                "First unmatched IDs:"
            )

            print(
                unmatched
                .head(20)
                .to_string(index=False)
            )


print()
print("=" * 100)

print(
    "QC complete."
)

print(
    "No lipid-lowering/statin exposure has been defined."
)

print(
    "No regression model has been fitted."
)

print()
