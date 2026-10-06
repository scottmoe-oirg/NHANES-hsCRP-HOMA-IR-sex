from __future__ import annotations

import argparse
from pathlib import Path
import sys
import tomllib

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data import build_combined_dataset
from src.survey import (
    survey_linear_regression,
    joint_wald_test,
    linear_contrast,
)


# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

def load_config(path: Path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def egfr_2021_ckd_epi(scr, age, is_male):
    """
    2021 CKD-EPI creatinine equation.

    scr:
        Serum creatinine in mg/dL.

    age:
        Age in years.

    is_male:
        1 = male
        0 = female
    """

    if pd.isna(scr) or pd.isna(age) or pd.isna(is_male):
        return np.nan

    if scr <= 0:
        return np.nan

    if is_male == 1:
        kappa = 0.9
        alpha = -0.302
        sex_factor = 1.0

    elif is_male == 0:
        kappa = 0.7
        alpha = -0.241
        sex_factor = 1.012

    else:
        return np.nan

    ratio = scr / kappa

    return (
        142.0
        * min(ratio, 1.0) ** alpha
        * max(ratio, 1.0) ** (-1.200)
        * (0.9938 ** age)
        * sex_factor
    )


def load_creatinine(
    renal_dir: Path,
):
    """
    Load serum creatinine from the two NHANES periods used
    in Study #2.
    """

    paths = [
        renal_dir / "BIOPRO_I.xpt",
        renal_dir / "P_BIOPRO.xpt",
    ]

    frames = []

    for path in paths:

        if not path.exists():
            raise FileNotFoundError(
                f"Missing renal file: {path}"
            )

        x = pd.read_sas(
            path,
            format="xport",
        )

        needed = [
            "SEQN",
            "LBXSCR",
        ]

        missing = [
            col
            for col in needed
            if col not in x.columns
        ]

        if missing:
            raise RuntimeError(
                f"{path.name} missing columns: {missing}"
            )

        frames.append(
            x[needed].copy()
        )

    renal = pd.concat(
        frames,
        ignore_index=True,
    )

    if renal["SEQN"].duplicated().any():
        raise RuntimeError(
            "Duplicate SEQN found in combined renal data."
        )

    return renal


def fit_model(
    df,
    domain,
    label,
    renal_term=None,
):

    base = (
        "log_hs_CRP ~ HOMA_IR * is_male "
        "+ C(period) * "
        "(LBXGH + waist + AGE + Non_HDL)"
    )

    if renal_term is None:
        formula = base
    else:
        formula = (
            base
            + f" + {renal_term}"
        )

    result = survey_linear_regression(
        df,
        formula,
        domain,
    )

    interaction_term = (
        "HOMA_IR:is_male"
    )

    test = joint_wald_test(
        result,
        [interaction_term],
    )

    female = linear_contrast(
        result,
        {
            "HOMA_IR": 1.0,
        },
    )

    male = linear_contrast(
        result,
        {
            "HOMA_IR": 1.0,
            interaction_term: 1.0,
        },
    )

    difference = linear_contrast(
        result,
        {
            interaction_term: 1.0,
        },
    )

    return {
        "model": label,
        "n": result["n"],
        "design_df": result["df"],

        "interaction_F": test["F"],
        "interaction_p": test["p"],

        "female_beta": female["beta"],
        "female_se": female["se"],
        "female_ci_low": female["ci_low"],
        "female_ci_high": female["ci_high"],
        "female_p": female["p"],

        "male_beta": male["beta"],
        "male_se": male["se"],
        "male_ci_low": male["ci_low"],
        "male_ci_high": male["ci_high"],
        "male_p": male["p"],

        "male_minus_female": difference["beta"],
        "difference_se": difference["se"],
        "difference_ci_low": difference["ci_low"],
        "difference_ci_high": difference["ci_high"],
        "difference_p": difference["p"],

        "formula": formula,
    }


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        default=str(
            ROOT
            / "config"
            / "analysis_config.toml"
        ),
    )

    parser.add_argument(
        "--data-dir",
        required=True,
    )

    parser.add_argument(
        "--renal-dir",
        default=str(
            ROOT
            / "data"
            / "raw_supplemental"
        ),
    )

    parser.add_argument(
        "--out-dir",
        default=str(
            ROOT
            / "output"
            / "sex_renal_sensitivity"
        ),
    )

    args = parser.parse_args()

    config = load_config(
        Path(args.config)
    )

    # --------------------------------------------------------
    # Build validated Study #2 dataset
    # --------------------------------------------------------

    df, domain, flow = build_combined_dataset(
        Path(args.data_dir),
        config,
    )

    n_locked = int(
        domain.sum()
    )

    if n_locked != 3079:
        raise RuntimeError(
            "Expected primary analytic domain "
            f"n=3079; got {n_locked}"
        )

    # --------------------------------------------------------
    # Load and merge serum creatinine
    # --------------------------------------------------------

    renal = load_creatinine(
        Path(args.renal_dir)
    )

    original_rows = len(df)

    df = df.merge(
        renal,
        on="SEQN",
        how="left",
        validate="one_to_one",
        sort=False,
    )

    if len(df) != original_rows:
        raise RuntimeError(
            "Renal merge changed dataset row count."
        )

    # --------------------------------------------------------
    # Construct frozen renal variables
    # --------------------------------------------------------

    df["eGFR_2021"] = df.apply(
        lambda row: egfr_2021_ckd_epi(
            row["LBXSCR"],
            row["AGE"],
            row["is_male"],
        ),
        axis=1,
    )

    df["log_LBXSCR"] = np.where(
        df["LBXSCR"] > 0,
        np.log(df["LBXSCR"]),
        np.nan,
    )

    # --------------------------------------------------------
    # Frozen QC checks
    # --------------------------------------------------------

    locked = df.loc[
        domain
    ].copy()

    n_creatinine = int(
        locked["LBXSCR"]
        .notna()
        .sum()
    )

    if n_creatinine != 3075:
        raise RuntimeError(
            "Expected serum creatinine available "
            f"for n=3075; got {n_creatinine}"
        )

    n_egfr = int(
        locked["eGFR_2021"]
        .notna()
        .sum()
    )

    n_log_creatinine = int(
        locked["log_LBXSCR"]
        .notna()
        .sum()
    )

    if n_egfr != 3075:
        raise RuntimeError(
            "Expected eGFR available for "
            f"n=3075; got {n_egfr}"
        )

    if n_log_creatinine != 3075:
        raise RuntimeError(
            "Expected log creatinine available for "
            f"n=3075; got {n_log_creatinine}"
        )

    # Confirm frozen distributional QC.
    raw_skew = locked["LBXSCR"].skew()
    log_skew = locked["log_LBXSCR"].skew()

    if not np.isclose(
        raw_skew,
        15.759526,
        atol=1e-5,
    ):
        raise RuntimeError(
            "Raw creatinine skewness does not reproduce "
            f"frozen QC: {raw_skew}"
        )

    if not np.isclose(
        log_skew,
        0.709910,
        atol=1e-5,
    ):
        raise RuntimeError(
            "Log creatinine skewness does not reproduce "
            f"frozen QC: {log_skew}"
        )

    # --------------------------------------------------------
    # Common renal domain
    #
    # ALL THREE models use exactly the same participants.
    # --------------------------------------------------------

    renal_complete = (
        df["eGFR_2021"].notna()
        & df["log_LBXSCR"].notna()
    )

    common_domain = (
        domain
        & renal_complete
    )

    n_common = int(
        common_domain.sum()
    )

    if n_common != 3075:
        raise RuntimeError(
            "Expected renal common domain "
            f"n=3075; got {n_common}"
        )

    # --------------------------------------------------------
    # Fit all three models on identical n=3075 sample
    # --------------------------------------------------------

    rows = []

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary model, common sample",
            renal_term=None,
        )
    )

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary + continuous eGFR",
            renal_term="eGFR_2021",
        )
    )

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary + log serum creatinine",
            renal_term="log_LBXSCR",
        )
    )

    results = pd.DataFrame(
        rows
    )

    # --------------------------------------------------------
    # Quantify change in focal sex-interaction estimate
    # --------------------------------------------------------

    reference = results.iloc[0]
    egfr_adjusted = results.iloc[1]
    logcr_adjusted = results.iloc[2]

    def change_metrics(adjusted):

        delta = (
            adjusted["male_minus_female"]
            - reference["male_minus_female"]
        )

        abs_ref = abs(
            reference["male_minus_female"]
        )

        abs_adjusted = abs(
            adjusted["male_minus_female"]
        )

        pct = (
            100
            * (abs_adjusted - abs_ref)
            / abs_ref
        )

        return delta, pct

    egfr_delta, egfr_pct = (
        change_metrics(
            egfr_adjusted
        )
    )

    logcr_delta, logcr_pct = (
        change_metrics(
            logcr_adjusted
        )
    )

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    out = Path(
        args.out_dir
    )

    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    results.to_csv(
        out / "renal_sensitivity.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Display
    # --------------------------------------------------------

    print()
    print(
        "STAGE 2D: RENAL FUNCTION SENSITIVITY"
    )
    print("=" * 115)

    print(
        "Locked primary analytic domain: "
        f"{n_locked:,}"
    )

    print(
        "Renal common domain: "
        f"{n_common:,}"
    )

    print(
        "Excluded for missing renal data: "
        f"{n_locked - n_common:,}"
    )

    print()
    print(
        "Frozen renal specifications:"
    )

    print(
        "  Principal:     continuous 2021 CKD-EPI "
        "creatinine eGFR"
    )

    print(
        "  Complementary: log serum creatinine "
        "(log LBXSCR)"
    )

    print()
    print(
        "Both renal specifications were frozen "
        "before outcome-model fitting."
    )

    print(
        "All three models below use exactly the same "
        f"n={n_common:,} participants."
    )

    print()
    print("RESULTS")
    print("-" * 115)

    display_cols = [
        "model",
        "n",
        "design_df",
        "interaction_F",
        "interaction_p",
        "female_beta",
        "male_beta",
        "male_minus_female",
        "difference_ci_low",
        "difference_ci_high",
    ]

    print(
        results[
            display_cols
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )

    print()
    print(
        "CHANGE IN MALE-MINUS-FEMALE CONTRAST"
    )
    print("-" * 115)

    print(
        "Common-sample reference:          "
        f"{reference['male_minus_female']:.6f}"
    )

    print()
    print(
        "After continuous eGFR adjustment: "
        f"{egfr_adjusted['male_minus_female']:.6f}"
    )

    print(
        "  Change in contrast:             "
        f"{egfr_delta:.6f}"
    )

    print(
        "  Percent change in absolute "
        "magnitude: "
        f"{egfr_pct:.2f}%"
    )

    print()
    print(
        "After log-creatinine adjustment:  "
        f"{logcr_adjusted['male_minus_female']:.6f}"
    )

    print(
        "  Change in contrast:             "
        f"{logcr_delta:.6f}"
    )

    print(
        "  Percent change in absolute "
        "magnitude: "
        f"{logcr_pct:.2f}%"
    )

    print()
    print("SEX-SPECIFIC SLOPES")
    print("-" * 115)

    for _, row in results.iterrows():

        print()
        print(
            row["model"]
        )

        print(
            "  Female beta = "
            f"{row['female_beta']:.6f}"
            "  "
            f"95% CI "
            f"[{row['female_ci_low']:.6f}, "
            f"{row['female_ci_high']:.6f}]"
        )

        print(
            "  Male beta   = "
            f"{row['male_beta']:.6f}"
            "  "
            f"95% CI "
            f"[{row['male_ci_low']:.6f}, "
            f"{row['male_ci_high']:.6f}]"
        )

    print()
    print("MODELS")
    print("-" * 115)

    for _, row in results.iterrows():

        print()
        print(
            row["model"]
        )

        print(
            row["formula"]
        )

    print()
    print(
        "No participants were excluded on the basis "
        "of low eGFR or high serum creatinine."
    )

    print(
        "eGFR thresholds were used for descriptive "
        "QC only."
    )

    print()


if __name__ == "__main__":
    main()
