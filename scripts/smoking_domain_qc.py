#!/usr/bin/env python3

from pathlib import Path
import sys
import tomllib

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data import build_combined_dataset


def load_config(path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def classify_smoking(row):

    ever = row["SMQ020"]
    current = row["SMQ040"]

    if ever == 2:
        return "Never"

    if ever == 1 and current == 3:
        return "Former"

    if ever == 1 and current in (1, 2):
        return "Current"

    return np.nan


config = load_config(
    ROOT / "config" / "analysis_config.toml"
)

data_dir = (
    ROOT.parent
    / "hsCRP-HOMA-IR"
    / "data_raw"
)

df, domain, flow = build_combined_dataset(
    data_dir,
    config,
)

if int(domain.sum()) != 3079:
    raise RuntimeError(
        f"Expected n=3079; got {int(domain.sum())}"
    )


# ------------------------------------------------------------
# Read supplemental smoking files
# ------------------------------------------------------------

supp = ROOT / "data" / "raw_supplemental"

files = [
    supp / "SMQ_I.xpt",
    supp / "P_SMQ.xpt",
]

frames = []

for path in files:

    smq = pd.read_sas(
        path,
        format="xport",
    )

    smq = smq[
        ["SEQN", "SMQ020", "SMQ040"]
    ].copy()

    frames.append(smq)


smoking = pd.concat(
    frames,
    ignore_index=True,
)

if smoking["SEQN"].duplicated().any():
    raise RuntimeError(
        "Duplicate SEQN found in smoking data."
    )

smoking["smoking_status"] = smoking.apply(
    classify_smoking,
    axis=1,
)


# ------------------------------------------------------------
# Merge smoking onto existing combined NHANES dataset
# ------------------------------------------------------------

x = df.merge(
    smoking[
        [
            "SEQN",
            "SMQ020",
            "SMQ040",
            "smoking_status",
        ]
    ],
    on="SEQN",
    how="left",
    validate="one_to_one",
)

# domain was created before merge; merge should preserve row order/count.
if len(x) != len(df):
    raise RuntimeError(
        "Merge changed number of rows."
    )


analytic = x.loc[domain].copy()

analytic["sex"] = analytic["is_male"].map(
    {
        0.0: "Female",
        1.0: "Male",
    }
)


print()
print("SMOKING QC WITHIN STUDY 2 ANALYTIC DOMAIN")
print("=" * 78)

print(
    f"\nLocked analytic domain: "
    f"{len(analytic):,}"
)


classified = (
    analytic["smoking_status"]
    .notna()
    .sum()
)

missing = (
    analytic["smoking_status"]
    .isna()
    .sum()
)

print(
    f"Classifiable smoking status: "
    f"{classified:,}"
)

print(
    f"Missing/unclassifiable:       "
    f"{missing:,}"
)

print(
    f"Percent classifiable:         "
    f"{100 * classified / len(analytic):.2f}%"
)


print("\nSMOKING STATUS — UNWEIGHTED")
print("-" * 78)

print(
    pd.crosstab(
        analytic["sex"],
        analytic["smoking_status"],
        margins=True,
        dropna=False,
    ).to_string()
)


print("\nSMOKING STATUS — ROW PERCENT WITHIN SEX")
print("-" * 78)

tab = pd.crosstab(
    analytic["sex"],
    analytic["smoking_status"],
    normalize="index",
)

print(
    (100 * tab)
    .round(2)
    .to_string()
)


print("\nMISSING / UNCLASSIFIABLE BY SEX")
print("-" * 78)

for sex in ["Female", "Male"]:

    z = analytic[
        analytic["sex"] == sex
    ]

    n_missing = (
        z["smoking_status"]
        .isna()
        .sum()
    )

    print(
        f"{sex:8s}: "
        f"{n_missing:,} / {len(z):,} "
        f"({100*n_missing/len(z):.2f}%)"
    )


print("\nRAW VALUES AMONG UNCLASSIFIABLE ANALYTIC PARTICIPANTS")
print("-" * 78)

bad = analytic[
    analytic["smoking_status"].isna()
][
    [
        "SEQN",
        "sex",
        "AGE",
        "SMQ020",
        "SMQ040",
    ]
]

if len(bad) == 0:
    print("None.")
else:
    print(
        bad.to_string(
            index=False
        )
    )


print()
print("No regression model has been fitted.")
print()
