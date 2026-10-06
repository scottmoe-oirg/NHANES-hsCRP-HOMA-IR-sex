#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np


ROOT = Path(__file__).resolve().parents[1]

SUPP = ROOT / "data" / "raw_supplemental"

FILES = {
    "2015-2016": SUPP / "SMQ_I.xpt",
    "2017-March2020": SUPP / "P_SMQ.xpt",
}


def classify_smoking(row):
    """
    Pre-specified smoking-status definition:

    Never:
        SMQ020 = 2

    Former:
        SMQ020 = 1 and SMQ040 = 3

    Current:
        SMQ020 = 1 and SMQ040 in {1, 2}

    Anything else:
        missing / unclassifiable
    """

    ever = row["SMQ020"]
    current = row["SMQ040"]

    if ever == 2:
        return "Never"

    if ever == 1:
        if current == 3:
            return "Former"
        if current in (1, 2):
            return "Current"

    return np.nan


frames = []

for period, path in FILES.items():

    if not path.exists():
        raise FileNotFoundError(
            f"Missing file: {path}"
        )

    df = pd.read_sas(
        path,
        format="xport",
    )

    needed = ["SEQN", "SMQ020", "SMQ040"]

    missing_columns = [
        x for x in needed
        if x not in df.columns
    ]

    if missing_columns:
        raise RuntimeError(
            f"{path.name} is missing columns: "
            f"{missing_columns}"
        )

    x = df[needed].copy()
    x["source_period"] = period

    x["smoking_status"] = x.apply(
        classify_smoking,
        axis=1,
    )

    frames.append(x)


smoking = pd.concat(
    frames,
    ignore_index=True,
)


print()
print("SMOKING QC")
print("=" * 72)

print("\nSOURCE FILES")
print("-" * 72)

for period, path in FILES.items():
    n = (
        smoking["source_period"]
        .eq(period)
        .sum()
    )

    print(
        f"{period:20s} "
        f"{path.name:12s} "
        f"n = {n:,}"
    )


print("\nDUPLICATE SEQN CHECK")
print("-" * 72)

duplicates = smoking["SEQN"].duplicated(
    keep=False
)

print(
    f"Duplicate SEQN across combined files: "
    f"{duplicates.sum():,}"
)


print("\nRAW SMQ020 COUNTS")
print("-" * 72)

print(
    pd.crosstab(
        smoking["source_period"],
        smoking["SMQ020"],
        dropna=False,
    ).to_string()
)


print("\nRAW SMQ040 COUNTS")
print("-" * 72)

print(
    pd.crosstab(
        smoking["source_period"],
        smoking["SMQ040"],
        dropna=False,
    ).to_string()
)


print("\nCONSTRUCTED SMOKING STATUS")
print("-" * 72)

status_table = pd.crosstab(
    smoking["source_period"],
    smoking["smoking_status"],
    dropna=False,
)

print(status_table.to_string())


print("\nOVERALL CONSTRUCTED STATUS")
print("-" * 72)

print(
    smoking["smoking_status"]
    .value_counts(
        dropna=False
    )
    .to_string()
)


print("\nQC OF STRUCTURAL SMQ040 MISSINGNESS")
print("-" * 72)

never = smoking["SMQ020"] == 2

print(
    "Participants with SMQ020=2 (never 100 cigarettes): "
    f"{never.sum():,}"
)

print(
    "Among them, SMQ040 missing: "
    f"{smoking.loc[never, 'SMQ040'].isna().sum():,}"
)


print("\nUNCLASSIFIABLE RESPONSES")
print("-" * 72)

unclassifiable = smoking[
    smoking["smoking_status"].isna()
]

print(
    f"Total unclassifiable: "
    f"{len(unclassifiable):,}"
)

if len(unclassifiable):
    print()
    print(
        unclassifiable[
            [
                "SEQN",
                "source_period",
                "SMQ020",
                "SMQ040",
            ]
        ]
        .head(30)
        .to_string(index=False)
    )


print()
print("No regression model has been fitted.")
print()
