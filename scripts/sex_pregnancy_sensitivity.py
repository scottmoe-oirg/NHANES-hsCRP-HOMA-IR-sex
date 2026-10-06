import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "raw"
SUPP = DATA_DIR
CONFIG_PATH = ROOT / "config" / "analysis_config.toml"
OUT_DIR = ROOT / "output" / "sex_pregnancy_sensitivity"

OUT_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT))

from run_analysis import load_config
from src.data import build_combined_dataset
from src.survey import survey_linear_regression, joint_wald_test, linear_contrast

# ------------------------------------------------------------
# Build the validated Study #2 domain
# ------------------------------------------------------------

config = load_config(CONFIG_PATH)
df, primary, audit = build_combined_dataset(DATA_DIR, config)

if int(primary.sum()) != 3079:
    raise RuntimeError(
        f"Expected locked Study #2 domain n=3079, "
        f"got {int(primary.sum())}."
    )

# ------------------------------------------------------------
# Add pregnancy status
# ------------------------------------------------------------

demo_i = pd.read_sas(SUPP / "DEMO_I.xpt", format="xport")
p_demo = pd.read_sas(SUPP / "P_DEMO.xpt", format="xport")

preg_i = demo_i[["SEQN", "RIDEXPRG"]].copy()
preg_p = p_demo[["SEQN", "RIDEXPRG"]].copy()

preg = pd.concat([preg_i, preg_p], ignore_index=True)

if preg["SEQN"].duplicated().any():
    raise RuntimeError("Duplicate SEQN in pregnancy data.")

# Merge pregnancy status into the full combined NHANES frame.
# This is important because the survey regression machinery uses
# domain score-zeroing rather than physically discarding all
# observations outside the analytic domain.

df = df.merge(
    preg,
    on="SEQN",
    how="left",
    validate="one_to_one",
)

if len(df) != 25531:
    raise RuntimeError(
        f"Pregnancy merge changed full-frame row count to {len(df)}."
    )

# Re-establish the primary domain mask after merge using SEQN.
primary_seqn = set(
    build_combined_dataset(DATA_DIR, config)[0]
    .loc[primary, "SEQN"]
)

domain_full = df["SEQN"].isin(primary_seqn)

if int(domain_full.sum()) != 3079:
    raise RuntimeError(
        f"Expected primary domain n=3079 after merge, "
        f"got {int(domain_full.sum())}."
    )

# ------------------------------------------------------------
# Frozen pregnancy sensitivity domains
# ------------------------------------------------------------

# Full locked primary domain.
domain_primary = domain_full.copy()

# Principal Stage 2F sensitivity:
# exclude confirmed current pregnancy only (RIDEXPRG == 1).
domain_no_pregnant = (
    domain_full
    & ~df["RIDEXPRG"].eq(1)
)

# Complementary conservative sensitivity:
# exclude confirmed pregnancy AND cannot-determine status.
#
# IMPORTANT:
# NaN RIDEXPRG values are retained. In the locked cohort these
# correspond to females age >=45 or participants for whom pregnancy
# status is not applicable; QC established that no females age
# 20-44 had NaN RIDEXPRG.
domain_no_preg_or_unknown = (
    domain_full
    & ~df["RIDEXPRG"].isin([1, 3])
)

# ------------------------------------------------------------
# QC domain sizes
# ------------------------------------------------------------

expected = {
    "full_primary": 3079,
    "exclude_pregnant": 3031,
    "exclude_pregnant_and_cannot_determine": 3016,
}

domains = {
    "full_primary": domain_primary,
    "exclude_pregnant": domain_no_pregnant,
    "exclude_pregnant_and_cannot_determine":
        domain_no_preg_or_unknown,
}

print("\n=== DOMAIN QC ===")

for name, mask in domains.items():
    n = int(mask.sum())

    female_n = int(
        (mask & df["is_male"].eq(0)).sum()
    )

    male_n = int(
        (mask & df["is_male"].eq(1)).sum()
    )

    print(
        f"{name}: n={n}, "
        f"female={female_n}, male={male_n}"
    )

    if n != expected[name]:
        raise RuntimeError(
            f"{name}: expected n={expected[name]}, got {n}"
        )

    if male_n != 1441:
        raise RuntimeError(
            f"{name}: expected male n=1441, got {male_n}"
        )

# ------------------------------------------------------------
# Locked Study #2 primary interaction model
# ------------------------------------------------------------

formula = (
    "log_hs_CRP ~ HOMA_IR * is_male "
    "+ C(period) * (LBXGH + waist + AGE + Non_HDL)"
)

# ------------------------------------------------------------
# Fit helper
# ------------------------------------------------------------

def fit_domain(label, domain_mask):

    result = survey_linear_regression(
        df,
        formula,
        domain_mask,
    )

    # Formal HOMA-IR x sex interaction.
    interaction = joint_wald_test(
        result,
        ["HOMA_IR:is_male"],
    )

    # Female slope:
    # female is reference group (is_male = 0).
    female = linear_contrast(
        result,
        {"HOMA_IR": 1.0},
    )

    # Male slope:
    # HOMA_IR + HOMA_IR:is_male.
    male = linear_contrast(
        result,
        {
            "HOMA_IR": 1.0,
            "HOMA_IR:is_male": 1.0,
        },
    )

    # Male-minus-female slope difference.
    difference = linear_contrast(
        result,
        {"HOMA_IR:is_male": 1.0},
    )

    row = {
        "analysis": label,
        "n": result["n"],
        "design_df": result["df"],

        "interaction_F": interaction["F"],
        "interaction_p": interaction["p"],

        "female_slope": female["beta"],
        "female_se": female["se"],
        "female_ci_low": female["ci_low"],
        "female_ci_high": female["ci_high"],
        "female_p": female["p"],

        "male_slope": male["beta"],
        "male_se": male["se"],
        "male_ci_low": male["ci_low"],
        "male_ci_high": male["ci_high"],
        "male_p": male["p"],

        "male_minus_female": difference["beta"],
        "difference_se": difference["se"],
        "difference_ci_low": difference["ci_low"],
        "difference_ci_high": difference["ci_high"],
        "difference_p": difference["p"],
    }

    return row

# ------------------------------------------------------------
# Run frozen analyses
# ------------------------------------------------------------

rows = []

for label, mask in domains.items():
    print(f"\nFitting: {label}")
    rows.append(
        fit_domain(label, mask)
    )

results = pd.DataFrame(rows)

# ------------------------------------------------------------
# Quantify change from full primary estimate
# ------------------------------------------------------------

primary_diff = results.loc[
    results["analysis"] == "full_primary",
    "male_minus_female",
].iloc[0]

results["change_from_full_primary"] = (
    results["male_minus_female"] - primary_diff
)

results["absolute_contrast_change_percent"] = (
    (
        results["male_minus_female"].abs()
        - abs(primary_diff)
    )
    / abs(primary_diff)
    * 100.0
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

out_file = OUT_DIR / "pregnancy_sensitivity.csv"
results.to_csv(out_file, index=False)

# ------------------------------------------------------------
# Print concise results
# ------------------------------------------------------------

print("\n=== STAGE 2F PREGNANCY SENSITIVITY ===")

for _, r in results.iterrows():

    print(f"\n{r['analysis']}")
    print(
        f"n={int(r['n'])}, "
        f"df={r['design_df']:.0f}"
    )

    print(
        f"Interaction: "
        f"F={r['interaction_F']:.6f}, "
        f"p={r['interaction_p']:.6f}"
    )

    print(
        "Female slope: "
        f"{r['female_slope']:.6f} "
        f"(95% CI "
        f"{r['female_ci_low']:.6f}, "
        f"{r['female_ci_high']:.6f})"
    )

    print(
        "Male slope: "
        f"{r['male_slope']:.6f} "
        f"(95% CI "
        f"{r['male_ci_low']:.6f}, "
        f"{r['male_ci_high']:.6f})"
    )

    print(
        "Male-minus-female: "
        f"{r['male_minus_female']:.6f} "
        f"(95% CI "
        f"{r['difference_ci_low']:.6f}, "
        f"{r['difference_ci_high']:.6f})"
    )

    print(
        "Change from full primary: "
        f"{r['change_from_full_primary']:.6f}"
    )

    print(
        "Change in absolute contrast magnitude: "
        f"{r['absolute_contrast_change_percent']:.2f}%"
    )

print("\nSaved:")
print(out_file)
