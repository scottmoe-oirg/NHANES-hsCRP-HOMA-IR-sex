#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys
import tomllib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data import build_combined_dataset

DATA_DIR = ROOT / "data" / "raw"

CONFIG = ROOT / "config" / "analysis_config.toml"

OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ------------------------------------------------------------
# Locked Study #2 cohort
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


locked["sex"] = np.where(
    locked["is_male"] == 1,
    "Male",
    "Female",
)

locked["high_hscrp"] = (
    locked["hs_CRP"] > 10.0
)


if int(locked["high_hscrp"].sum()) != 184:
    raise RuntimeError(
        "Expected 184 participants with hsCRP >10 mg/L."
    )


# ------------------------------------------------------------
# Basic QC
# ------------------------------------------------------------

print()
print("SCATTERPLOT QC")
print("=" * 70)

for sex in ["Female", "Male"]:

    g = locked.loc[
        locked["sex"] == sex
    ]

    print(
        f"{sex:<8} "
        f"n={len(g):4d}  "
        f"hsCRP>10={int(g['high_hscrp'].sum()):3d}"
    )


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig, axes = plt.subplots(
    1,
    2,
    figsize=(13, 6),
    sharex=True,
    sharey=True,
)


for ax, sex in zip(
    axes,
    ["Female", "Male"],
):

    g = locked.loc[
        locked["sex"] == sex
    ]

    low = g.loc[
        ~g["high_hscrp"]
    ]

    high = g.loc[
        g["high_hscrp"]
    ]

    # hsCRP <= 10
    ax.scatter(
        low["HOMA_IR"],
        low["log_hs_CRP"],
        s=16,
        alpha=0.25,
        label="hsCRP ≤10 mg/L",
    )

    # hsCRP > 10
    ax.scatter(
        high["HOMA_IR"],
        high["log_hs_CRP"],
        s=38,
        alpha=0.80,
        marker="o",
        edgecolors="black",
        linewidths=0.5,
        label="hsCRP >10 mg/L",
    )

    # Diagnostic display zoom only.
    # No observations are removed from the analytic dataset.
    ax.set_xlim(0, 12)

    # Threshold corresponding to hsCRP = 10
    ax.axhline(
        np.log(10),
        linestyle="--",
        linewidth=1.2,
        alpha=0.8,
    )

    ax.set_title(
        f"{sex} (n={len(g):,})"
    )

    ax.set_xlabel(
        "HOMA-IR"
    )

    ax.grid(
        alpha=0.15
    )


axes[0].set_ylabel(
    "log(hsCRP)"
)


axes[1].legend(
    frameon=False,
    loc="upper right",
)


fig.suptitle(
    "HOMA-IR and hsCRP by Sex: Zoomed Diagnostic View",
    fontsize=15,
)

fig.text(
    0.5,
    0.01,
    "Display restricted to HOMA-IR 0–12 for visualization only; "
    "no observations were excluded from analysis. "
    "Dashed line corresponds to hsCRP = 10 mg/L.",
    ha="center",
    fontsize=9,
)

plt.tight_layout(
    rect=[0, 0.05, 1, 0.95]
)

png_path = (
    OUT_DIR
    / "homa_hscrp_scatter_by_sex_zoom.png"
)

pdf_path = (
    OUT_DIR
    / "homa_hscrp_scatter_by_sex_zoom.pdf"
)

fig.savefig(
    png_path,
    dpi=300,
    bbox_inches="tight",
)

fig.savefig(
    pdf_path,
    bbox_inches="tight",
)


print()
print("SAVED")
print("-" * 70)
print(png_path)
print(pdf_path)

print()
print(
    "No fitted lines, hypothesis tests, or new models were used."
)


plt.show()
