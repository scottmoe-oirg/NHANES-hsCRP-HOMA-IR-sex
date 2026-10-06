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
SUPP_DIR = DATA_DIR

CONFIG = ROOT / "config" / "analysis_config.toml"


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def weighted_quantile(values, weights, q):
    """
    Simple weighted empirical quantile.
    Descriptive only; no variance estimation.
    """

    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)

    keep = (
        np.isfinite(values)
        & np.isfinite(weights)
        & (weights > 0)
    )

    values = values[keep]
    weights = weights[keep]

    if len(values) == 0:
        return np.nan

    order = np.argsort(values)

    values = values[order]
    weights = weights[order]

    cumulative = np.cumsum(weights)

    cutoff = q * cumulative[-1]

    idx = np.searchsorted(
        cumulative,
        cutoff,
        side="left",
    )

    idx = min(idx, len(values) - 1)

    return float(values[idx])


def describe_unweighted(x):
    x = pd.to_numeric(
        x,
        errors="coerce",
    ).dropna()

    return {
        "n": len(x),
        "mean": x.mean(),
        "median": x.median(),
        "q25": x.quantile(0.25),
        "q75": x.quantile(0.75),
        "p90": x.quantile(0.90),
    }


def describe_weighted(x, w):
    tmp = pd.DataFrame(
        {
            "x": pd.to_numeric(
                x,
                errors="coerce",
            ),
            "w": pd.to_numeric(
                w,
                errors="coerce",
            ),
        }
    )

    tmp = tmp.loc[
        tmp["x"].notna()
        & tmp["w"].notna()
        & (tmp["w"] > 0)
    ]

    if len(tmp) == 0:
        return {
            "weighted_mean": np.nan,
            "weighted_median": np.nan,
            "weighted_q25": np.nan,
            "weighted_q75": np.nan,
            "weighted_p90": np.nan,
        }

    weighted_mean = np.average(
        tmp["x"],
        weights=tmp["w"],
    )

    return {
        "weighted_mean": weighted_mean,
        "weighted_median": weighted_quantile(
            tmp["x"],
            tmp["w"],
            0.50,
        ),
        "weighted_q25": weighted_quantile(
            tmp["x"],
            tmp["w"],
            0.25,
        ),
        "weighted_q75": weighted_quantile(
            tmp["x"],
            tmp["w"],
            0.75,
        ),
        "weighted_p90": weighted_quantile(
            tmp["x"],
            tmp["w"],
            0.90,
        ),
    }


# ------------------------------------------------------------
# Build locked Study #2 cohort
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
        f"Expected locked domain n=3079; "
        f"got {len(locked)}"
    )


# ------------------------------------------------------------
# Identify survey weight used by validated dataset
# ------------------------------------------------------------

weight_candidates = [
    c for c in locked.columns
    if (
        "WTSAF" in c.upper()
        or c.lower() == "weight"
        or "fast" in c.lower()
        and "weight" in c.lower()
    )
]

print()
print("POTENTIAL SURVEY WEIGHT COLUMNS")
print("-" * 100)

for c in weight_candidates:
    print(c)


# We expect WTSAF2YR to be the fasting examination weight
# underlying the combined-period analysis. If the validated
# dataset uses another derived weight column, stop rather than
# silently choosing the wrong one.

if "WTSAF2YR" in locked.columns:
    weight_col = "WTSAF2YR"
else:
    raise RuntimeError(
        "WTSAF2YR not found. Inspect candidate weight columns "
        "before proceeding."
    )


# ------------------------------------------------------------
# Read and combine CBC
# ------------------------------------------------------------

cbc_1516 = pd.read_sas(
    SUPP_DIR / "CBC_I.xpt",
    format="xport",
)

cbc_1720 = pd.read_sas(
    SUPP_DIR / "P_CBC.xpt",
    format="xport",
)

cbc = pd.concat(
    [cbc_1516, cbc_1720],
    ignore_index=True,
    sort=False,
)

if cbc["SEQN"].duplicated().any():
    raise RuntimeError(
        "Duplicate CBC SEQN values detected."
    )


merged = locked.merge(
    cbc,
    on="SEQN",
    how="left",
    validate="one_to_one",
)


# ------------------------------------------------------------
# Frozen hsCRP grouping
# ------------------------------------------------------------

merged["high_hscrp"] = (
    merged["hs_CRP"] > 10.0
)

if int(merged["high_hscrp"].sum()) != 184:
    raise RuntimeError(
        "Expected exactly 184 participants with hsCRP >10."
    )


merged["sex"] = np.where(
    merged["is_male"] == 1,
    "Male",
    "Female",
)

merged["hscrp_group"] = np.where(
    merged["high_hscrp"],
    ">10",
    "<=10",
)


# ------------------------------------------------------------
# Frozen CBC panel
# ------------------------------------------------------------

required = [
    "LBXWBCSI",   # total WBC
    "LBDNENO",    # absolute neutrophils
    "LBXNEPCT",   # neutrophil %
    "LBDLYMNO",   # absolute lymphocytes
    "LBDMONO",    # absolute monocytes
]

for c in required:
    if c not in merged.columns:
        raise RuntimeError(
            f"Required CBC variable missing: {c}"
        )


# NLR is descriptive only.

merged["NLR"] = np.where(
    merged["LBDLYMNO"] > 0,
    merged["LBDNENO"]
    / merged["LBDLYMNO"],
    np.nan,
)


panel = {
    "hsCRP": "hs_CRP",
    "HOMA_IR": "HOMA_IR",
    "WBC": "LBXWBCSI",
    "Absolute_neutrophils": "LBDNENO",
    "Neutrophil_percent": "LBXNEPCT",
    "Absolute_lymphocytes": "LBDLYMNO",
    "Absolute_monocytes": "LBDMONO",
    "NLR": "NLR",
}


# ------------------------------------------------------------
# Sex-specific frequency of hsCRP >10
# ------------------------------------------------------------

print()
print("HIGH-hsCRP CBC DESCRIPTIVE DIAGNOSTIC")
print("=" * 100)

print()
print("FREQUENCY OF hsCRP >10 mg/L BY SEX")
print("-" * 100)

for sex in ["Female", "Male"]:

    g = merged.loc[
        merged["sex"] == sex
    ]

    high = g["high_hscrp"]

    n_total = len(g)
    n_high = int(high.sum())

    unweighted_pct = (
        100.0 * n_high / n_total
    )

    w = g[weight_col]

    weighted_pct = (
        100.0
        * np.sum(
            w * high.astype(float)
        )
        / np.sum(w)
    )

    print(
        f"{sex:<8} "
        f"n={n_total:4d}  "
        f">10={n_high:3d}  "
        f"unweighted={unweighted_pct:6.2f}%  "
        f"weighted={weighted_pct:6.2f}%"
    )


# ------------------------------------------------------------
# Four-cell descriptive summaries
# ------------------------------------------------------------

rows = []

for sex in ["Female", "Male"]:

    for hgroup in ["<=10", ">10"]:

        g = merged.loc[
            (merged["sex"] == sex)
            & (merged["hscrp_group"] == hgroup)
        ]

        for label, variable in panel.items():

            unweighted = describe_unweighted(
                g[variable]
            )

            weighted = describe_weighted(
                g[variable],
                g[weight_col],
            )

            rows.append(
                {
                    "sex": sex,
                    "hscrp_group": hgroup,
                    "variable": label,
                    **unweighted,
                    **weighted,
                }
            )


summary = pd.DataFrame(rows)


print()
print("FOUR-CELL DESCRIPTIVE SUMMARIES")
print("=" * 100)

for sex in ["Female", "Male"]:

    for hgroup in ["<=10", ">10"]:

        print()
        print(
            f"{sex} | hsCRP {hgroup} mg/L"
        )

        print("-" * 100)

        block = summary.loc[
            (summary["sex"] == sex)
            & (
                summary["hscrp_group"]
                == hgroup
            )
        ]

        cols = [
            "variable",
            "n",
            "mean",
            "median",
            "q25",
            "q75",
            "p90",
            "weighted_mean",
            "weighted_median",
        ]

        print(
            block[cols].to_string(
                index=False,
                float_format=lambda x: f"{x:.3f}",
            )
        )


# ------------------------------------------------------------
# Focused look at the 184 only
# ------------------------------------------------------------

print()
print("FOCUSED SUMMARY: hsCRP >10 mg/L ONLY")
print("=" * 100)

high = merged.loc[
    merged["high_hscrp"]
].copy()

for sex in ["Female", "Male"]:

    g = high.loc[
        high["sex"] == sex
    ]

    print()
    print(
        f"{sex}: n={len(g)}"
    )

    print("-" * 100)

    for label, variable in panel.items():

        x = g[variable]

        print(
            f"{label:<24} "
            f"median={x.median():8.3f}  "
            f"IQR=[{x.quantile(.25):8.3f}, "
            f"{x.quantile(.75):8.3f}]  "
            f"max={x.max():8.3f}"
        )


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

out_dir = (
    ROOT
    / "output"
    / "hscrp_high_cbc_descriptive"
)

out_dir.mkdir(
    parents=True,
    exist_ok=True,
)

summary.to_csv(
    out_dir
    / "cbc_four_cell_summary.csv",
    index=False,
)

high[
    [
        "SEQN",
        "sex",
        "hs_CRP",
        "HOMA_IR",
        "LBXWBCSI",
        "LBDNENO",
        "LBXNEPCT",
        "LBDLYMNO",
        "LBDMONO",
        "NLR",
    ]
].to_csv(
    out_dir
    / "high_hscrp_individual_values.csv",
    index=False,
)


print()
print("OUTPUT")
print("-" * 100)

print(
    out_dir
    / "cbc_four_cell_summary.csv"
)

print(
    out_dir
    / "high_hscrp_individual_values.csv"
)

print()
print(
    "DESCRIPTIVE DIAGNOSTIC ONLY."
)

print(
    "No hypothesis tests, infection classifications, "
    "new exclusions, or HOMA-IR x sex models were fitted."
)
