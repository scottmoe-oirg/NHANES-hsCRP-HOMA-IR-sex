#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys
import tomllib

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data import build_combined_dataset

DATA_DIR = ROOT / "data" / "raw"
SUPP_DIR = DATA_DIR

CONFIG = ROOT / "config" / "analysis_config.toml"


def read_xpt(path):
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_sas(
        path,
        format="xport",
    )

    return df


# ------------------------------------------------------------
# Build locked Study #2 analytic cohort
# ------------------------------------------------------------

with open(CONFIG, "rb") as f:
    config = tomllib.load(f)


df, domain, flow = build_combined_dataset(
    DATA_DIR,
    config,
)

locked = df.loc[domain].copy()

if len(locked) != 3079:
    raise RuntimeError(
        f"Expected locked domain n=3079; got {len(locked)}"
    )


# ------------------------------------------------------------
# Read CBC files
# ------------------------------------------------------------

cbc_1516 = read_xpt(
    SUPP_DIR / "CBC_I.xpt"
)

cbc_1720 = read_xpt(
    SUPP_DIR / "P_CBC.xpt"
)


print()
print("HIGH-hsCRP CBC QC")
print("=" * 90)

print()
print("SOURCE FILES")
print("-" * 90)

print(
    f"CBC_I.xpt:   {len(cbc_1516):,} rows, "
    f"{len(cbc_1516.columns):,} variables"
)

print(
    f"P_CBC.xpt:   {len(cbc_1720):,} rows, "
    f"{len(cbc_1720.columns):,} variables"
)


# ------------------------------------------------------------
# Inventory variables
# ------------------------------------------------------------

print()
print("2015-2016 CBC VARIABLES")
print("-" * 90)

for c in cbc_1516.columns:
    print(c)


print()
print("2017-MARCH 2020 CBC VARIABLES")
print("-" * 90)

for c in cbc_1720.columns:
    print(c)


# ------------------------------------------------------------
# Common variables
# ------------------------------------------------------------

common = sorted(
    set(cbc_1516.columns)
    & set(cbc_1720.columns)
)

print()
print("COMMON CBC VARIABLES")
print("-" * 90)

for c in common:
    print(c)


# ------------------------------------------------------------
# Combine cycles and merge to locked cohort
# ------------------------------------------------------------

cbc = pd.concat(
    [cbc_1516, cbc_1720],
    ignore_index=True,
    sort=False,
)

if cbc["SEQN"].duplicated().any():
    dup = int(
        cbc["SEQN"].duplicated().sum()
    )

    raise RuntimeError(
        f"Duplicate CBC SEQN values detected: {dup}"
    )


merged = locked.merge(
    cbc,
    on="SEQN",
    how="left",
    validate="one_to_one",
)


# ------------------------------------------------------------
# Define frozen high-hsCRP group
# ------------------------------------------------------------

merged["high_hscrp"] = (
    merged["hs_CRP"] > 10.0
)

n_high = int(
    merged["high_hscrp"].sum()
)

if n_high != 184:
    raise RuntimeError(
        f"Expected 184 participants with hsCRP >10; "
        f"got {n_high}"
    )


print()
print("LOCKED COHORT")
print("-" * 90)

print(
    f"Total:                 {len(merged):,}"
)

print(
    f"hsCRP <=10 mg/L:       "
    f"{(~merged['high_hscrp']).sum():,}"
)

print(
    f"hsCRP >10 mg/L:        "
    f"{merged['high_hscrp'].sum():,}"
)


# ------------------------------------------------------------
# Sex distribution of the 184
# ------------------------------------------------------------

print()
print("SEX DISTRIBUTION: hsCRP >10 mg/L")
print("-" * 90)

sex_counts = (
    merged.loc[
        merged["high_hscrp"],
        "is_male",
    ]
    .value_counts(
        dropna=False
    )
    .sort_index()
)

print(
    "Female (is_male=0):",
    int(sex_counts.get(0.0, 0)),
)

print(
    "Male   (is_male=1):",
    int(sex_counts.get(1.0, 0)),
)


# ------------------------------------------------------------
# Coverage for common CBC variables among the 184
# ------------------------------------------------------------

high = merged.loc[
    merged["high_hscrp"]
].copy()


print()
print("CBC COVERAGE AMONG THE 184")
print("-" * 90)

for c in common:

    if c == "SEQN":
        continue

    n = int(
        high[c].notna().sum()
    )

    if n > 0:
        print(
            f"{c:<15} {n:>4}/184 "
            f"({100*n/184:6.2f}%)"
        )


# ------------------------------------------------------------
# Identify likely immune variables by name only.
# No biological thresholds or outcome modeling here.
# ------------------------------------------------------------

keywords = (
    "WBC",
    "NEU",
    "LYM",
    "MON",
    "EOS",
    "BAS",
)

candidate_vars = [
    c for c in common
    if any(
        key in c.upper()
        for key in keywords
    )
]


print()
print("CANDIDATE WHITE-CELL / DIFFERENTIAL VARIABLES")
print("-" * 90)

for c in candidate_vars:
    print(c)


print()
print("NO INFERENTIAL MODEL WAS FIT.")
print(
    "This script inventories CBC availability, linkage, "
    "sex composition, and coverage among the frozen "
    "hsCRP >10 mg/L group only."
)
