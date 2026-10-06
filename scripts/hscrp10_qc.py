#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys
import tomllib

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data import build_combined_dataset

DATA_DIR = ROOT / "data" / "raw"

CONFIG = ROOT / "config" / "analysis_config.toml"


with open(CONFIG, "rb") as f:
    config = tomllib.load(f)


df, domain, flow = build_combined_dataset(
    DATA_DIR,
    config,
)

if int(domain.sum()) != 3079:
    raise RuntimeError(
        f"Expected locked primary domain n=3079; "
        f"got {int(domain.sum())}"
    )


x = df.loc[domain].copy()


# ------------------------------------------------------------
# Locate the raw hsCRP variable used to construct log_hs_CRP
# ------------------------------------------------------------

print()
print("STAGE 2H PRE-OUTCOME QC: hsCRP <= 10 mg/L")
print("=" * 80)

print(
    f"Locked primary analytic domain: {len(x):,}"
)

print()
print("Potential hsCRP-related columns:")
print("-" * 80)

hscrp_cols = [
    c for c in x.columns
    if (
        "CRP" in c.upper()
        or "HSCRP" in c.upper()
    )
]

for c in hscrp_cols:
    print(c)


print()
print("Basic values for candidate columns:")
print("-" * 80)

for c in hscrp_cols:

    if pd.api.types.is_numeric_dtype(x[c]):

        z = x[c]

        print()
        print(c)

        print(
            f"  nonmissing = {z.notna().sum():,}"
        )

        print(
            f"  <= 0       = "
            f"{((z <= 0) & z.notna()).sum():,}"
        )

        print(
            f"  > 10       = "
            f"{((z > 10) & z.notna()).sum():,}"
        )

        if z.notna().any():

            print(
                f"  min        = {z.min():.6f}"
            )

            print(
                f"  median     = {z.median():.6f}"
            )

            print(
                f"  max        = {z.max():.6f}"
            )


print()
print("No HOMA-IR x sex outcome model was fitted.")
print(
    "This script performs hsCRP restriction QC only."
)
