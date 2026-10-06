from pathlib import Path
import sys

import numpy as np
import pandas as pd

# Allow imports from project root when script is run from scripts/
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.data import build_combined_dataset
from run_analysis import load_config


DATA_DIR = ROOT / "data" / "raw"
CONFIG = ROOT / "config" / "analysis_config.toml"
OUT_DIR = ROOT / "output" / "table1"

OUT_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Weighted descriptive functions
# ------------------------------------------------------------

def weighted_mean_se(x, w, strata, psu, domain):
    """
    Survey-weighted domain mean and Taylor-linearized SE.

    The mean is estimated for observations inside the domain, while the
    variance calculation retains the full survey design by assigning zero
    linearized contributions to observations outside the domain.
    """

    valid = (
        x.notna()
        & w.notna()
        & (w > 0)
        & strata.notna()
        & psu.notna()
    )

    d = domain.astype(bool) & valid

    if d.sum() == 0:
        return np.nan, np.nan

    # Domain weighted mean
    numerator = np.sum(w.loc[d] * x.loc[d])
    denominator = np.sum(w.loc[d])
    mu = numerator / denominator

    # Linearized contribution for a ratio/domain mean:
    #
    #     z_i = I(domain_i) * w_i * (x_i - mu) / sum_domain(w)
    #
    # Observations outside the domain remain in the survey design but
    # contribute zero to the estimating equation.
    z = pd.Series(0.0, index=x.index)

    z.loc[d] = (
        w.loc[d]
        * (x.loc[d] - mu)
        / denominator
    )

    design = valid | (~domain.astype(bool))

    tmp = pd.DataFrame({
        "stratum": strata.loc[design],
        "psu": psu.loc[design],
        "z": z.loc[design],
    })

    # Remove records lacking usable survey-design identifiers.
    tmp = tmp.dropna(subset=["stratum", "psu"])

    # Sum linearized contributions within each PSU.
    psu_totals = (
        tmp.groupby(
            ["stratum", "psu"],
            observed=True
        )["z"]
        .sum()
        .reset_index()
    )

    variance = 0.0

    for _, g in psu_totals.groupby(
        "stratum",
        observed=True
    ):
        m_h = len(g)

        if m_h <= 1:
            continue

        vals = g["z"].to_numpy()
        mean_h = vals.mean()

        variance += (
            m_h / (m_h - 1.0)
        ) * np.sum(
            (vals - mean_h) ** 2
        )

    se = np.sqrt(variance)

    return mu, se

def weighted_quantile(x, w, q):
    mask = x.notna() & w.notna() & (w > 0)

    x = x.loc[mask].astype(float).to_numpy()
    w = w.loc[mask].astype(float).to_numpy()

    order = np.argsort(x)
    x = x[order]
    w = w[order]

    cum_w = np.cumsum(w)
    cutoff = q * np.sum(w)

    idx = np.searchsorted(cum_w, cutoff, side="left")
    idx = min(idx, len(x) - 1)

    return x[idx]


def weighted_median_iqr(x, w):
    q25 = weighted_quantile(x, w, 0.25)
    q50 = weighted_quantile(x, w, 0.50)
    q75 = weighted_quantile(x, w, 0.75)

    return q50, q25, q75


# ------------------------------------------------------------
# Formatting
# ------------------------------------------------------------

def fmt_mean_se(mean, se, decimals=1):
    return f"{mean:.{decimals}f} ({se:.{decimals}f})"


def fmt_median_iqr(median, q25, q75, decimals=2):
    return (
        f"{median:.{decimals}f} "
        f"({q25:.{decimals}f}-{q75:.{decimals}f})"
    )


# ------------------------------------------------------------
# Build locked Study #2 cohort
# ------------------------------------------------------------

config = load_config(CONFIG)

df, domain, flow = build_combined_dataset(DATA_DIR, config)

analytic = df.loc[domain].copy()

print("=" * 72)
print("TABLE 1: LOCKED COHORT CHECK")
print("=" * 72)
print(f"Analytic n = {len(analytic):,}")

if len(analytic) != 3079:
    raise RuntimeError(
        f"Expected locked analytic n=3079, got {len(analytic)}"
    )

print("\nSex counts:")
print(
    analytic["is_male"]
    .value_counts(dropna=False)
    .sort_index()
)

if analytic["is_male"].isna().any():
    raise RuntimeError("Missing sex values found in locked analytic cohort.")


# ------------------------------------------------------------
# Population estimates
# ------------------------------------------------------------

sex_labels = {
    0.0: "Female",
    1.0: "Male",
}

population_rows = []

for sex_code, sex_label in sex_labels.items():
    g = analytic.loc[analytic["is_male"] == sex_code]

    population_rows.append({
        "sex": sex_label,
        "unweighted_n": len(g),
        "weighted_population": g["analysis_weight"].sum(),
        "weighted_population_millions":
            g["analysis_weight"].sum() / 1_000_000,
    })

population = pd.DataFrame(population_rows)

print("\nWeighted population representation:")
print(population.to_string(index=False))


# ------------------------------------------------------------
# Table 1 variable definitions
# ------------------------------------------------------------

mean_variables = [
    ("Age, years", "AGE", 1),
    ("HbA1c, %", "LBXGH", 2),
    ("Waist circumference, cm", "waist", 1),
    ("Non-HDL cholesterol, mg/dL", "Non_HDL", 1),
    ("Fasting glucose, mg/dL", "LBXGLU", 1),
]

median_variables = [
    (r"Fasting insulin, $\mu$U/mL", "LBXIN", 2),
    ("HOMA-IR", "HOMA_IR", 2),
    ("hsCRP, mg/L", "hs_CRP", 2),
]


# ------------------------------------------------------------
# Calculate descriptive statistics
# ------------------------------------------------------------

rows = []

for label, var, decimals in mean_variables:
    row = {
        "Characteristic": label,
        "Statistic": "Weighted mean (SE)",
    }

    for sex_code, sex_label in sex_labels.items():
        g = analytic.loc[analytic["is_male"] == sex_code]

        sex_domain = analytic["is_male"].eq(sex_code)

        mean, se = weighted_mean_se(
            analytic[var],
            analytic["analysis_weight"],
            analytic["stratum_u"],
            analytic["psu_u"],
            sex_domain,
)

        row[sex_label] = fmt_mean_se(mean, se, decimals)

    rows.append(row)


for label, var, decimals in median_variables:
    row = {
        "Characteristic": label,
        "Statistic": "Weighted median (IQR)",
    }

    for sex_code, sex_label in sex_labels.items():
        g = analytic.loc[analytic["is_male"] == sex_code]

        med, q25, q75 = weighted_median_iqr(
            g[var],
            g["analysis_weight"],
        )

        row[sex_label] = fmt_median_iqr(
            med, q25, q75, decimals
        )

    rows.append(row)


table1 = pd.DataFrame(rows)


# ------------------------------------------------------------
# QC output
# ------------------------------------------------------------

print("\n" + "=" * 72)
print("TABLE 1 DESCRIPTIVE STATISTICS")
print("=" * 72)
print(table1.to_string(index=False))


# ------------------------------------------------------------
# Save machine-readable outputs
# ------------------------------------------------------------

population.to_csv(
    OUT_DIR / "table1_population_by_sex.csv",
    index=False,
)

table1.to_csv(
    OUT_DIR / "table1_descriptive_statistics.csv",
    index=False,
)


# ------------------------------------------------------------
# Create LaTeX table
# ------------------------------------------------------------

female_n = int(
    population.loc[
        population["sex"] == "Female",
        "unweighted_n"
    ].iloc[0]
)

male_n = int(
    population.loc[
        population["sex"] == "Male",
        "unweighted_n"
    ].iloc[0]
)

female_pop = float(
    population.loc[
        population["sex"] == "Female",
        "weighted_population_millions"
    ].iloc[0]
)

male_pop = float(
    population.loc[
        population["sex"] == "Male",
        "weighted_population_millions"
    ].iloc[0]
)


latex_lines = []

latex_lines.append(r"\begin{table}[htbp]")
latex_lines.append(r"\centering")
latex_lines.append(
    r"\caption{Characteristics of the primary analytic population by sex}"
)
latex_lines.append(r"\label{tab:table1}")
latex_lines.append(r"\begin{threeparttable}")
latex_lines.append(r"\begin{tabular}{lcc}")
latex_lines.append(r"\toprule")
latex_lines.append(
    f"Characteristic & Female ($n={female_n:,}$) "
    f"& Male ($n={male_n:,}$) \\\\"
)
latex_lines.append(r"\midrule")

latex_lines.append(
    f"Weighted population, millions "
    f"& {female_pop:.1f} & {male_pop:.1f} \\\\"
)

for _, r in table1.iterrows():
    label = r["Characteristic"]
    female = r["Female"]
    male = r["Male"]

    # Escape percent sign for LaTeX.
    label = label.replace("%", r"\%")

    latex_lines.append(
        f"{label} & {female} & {male} \\\\"
    )

latex_lines.append(r"\bottomrule")
latex_lines.append(r"\end{tabular}")

latex_lines.append(r"\begin{tablenotes}[flushleft]")
latex_lines.append(r"\footnotesize")
latex_lines.append(
    r"\item Values are survey-weighted mean (SE) for age, "
    r"HbA1c, waist circumference, non-HDL cholesterol, and "
    r"fasting glucose, and survey-weighted median (IQR) for "
    r"fasting insulin, HOMA-IR, and hsCRP."
)
latex_lines.append(
    r"\item The analytic population included adults aged "
    r"$\geq20$ years without diagnosed diabetes and with "
    r"HbA1c $<5.7\%$, meeting the fasting-subsample and "
    r"complete-case criteria for the primary model."
)
latex_lines.append(
    r"\item Weighted population estimates use the combined "
    r"NHANES 2015--2016 and 2017--March 2020 pre-pandemic "
    r"fasting-subsample weights."
)
latex_lines.append(
    r"\item No hypothesis tests comparing baseline characteristics "
    r"by sex were performed."
)
latex_lines.append(r"\end{tablenotes}")
latex_lines.append(r"\end{threeparttable}")
latex_lines.append(r"\end{table}")

latex_text = "\n".join(latex_lines)

latex_path = OUT_DIR / "table1_by_sex.tex"
latex_path.write_text(latex_text)

print("\n" + "=" * 72)
print("OUTPUT FILES")
print("=" * 72)
print(OUT_DIR / "table1_population_by_sex.csv")
print(OUT_DIR / "table1_descriptive_statistics.csv")
print(latex_path)

print("\nLaTeX preview:\n")
print(latex_text)
