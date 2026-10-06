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


def load_alt(
    hepatic_dir: Path,
):
    """
    Load ALT (LBXSATSI) from the two NHANES periods used
    in Study #2.
    """

    paths = [
        hepatic_dir / "BIOPRO_I.xpt",
        hepatic_dir / "P_BIOPRO.xpt",
    ]

    frames = []

    for path in paths:

        if not path.exists():
            raise FileNotFoundError(
                f"Missing hepatic file: {path}"
            )

        x = pd.read_sas(
            path,
            format="xport",
        )

        needed = [
            "SEQN",
            "LBXSATSI",
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

    hepatic = pd.concat(
        frames,
        ignore_index=True,
    )

    if hepatic["SEQN"].duplicated().any():
        raise RuntimeError(
            "Duplicate SEQN found in combined hepatic data."
        )

    return hepatic


def fit_model(
    df,
    domain,
    label,
    add_alt=False,
):

    base = (
        "log_hs_CRP ~ HOMA_IR * is_male "
        "+ C(period) * "
        "(LBXGH + waist + AGE + Non_HDL)"
    )

    if add_alt:
        formula = (
            base
            + " + log_ALT"
        )
    else:
        formula = base

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
        "--hepatic-dir",
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
            / "sex_hepatic_sensitivity"
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
    # Load and merge ALT
    # --------------------------------------------------------

    hepatic = load_alt(
        Path(args.hepatic_dir)
    )

    original_rows = len(df)

    df = df.merge(
        hepatic,
        on="SEQN",
        how="left",
        validate="one_to_one",
        sort=False,
    )

    if len(df) != original_rows:
        raise RuntimeError(
            "ALT merge changed dataset row count."
        )

    # --------------------------------------------------------
    # Construct frozen log-ALT variable
    # --------------------------------------------------------

    df["log_ALT"] = np.where(
        df["LBXSATSI"] > 0,
        np.log(df["LBXSATSI"]),
        np.nan,
    )

    # --------------------------------------------------------
    # Reproduce frozen QC
    # --------------------------------------------------------

    locked = df.loc[
        domain
    ].copy()

    n_alt = int(
        locked["LBXSATSI"]
        .notna()
        .sum()
    )

    n_log_alt = int(
        locked["log_ALT"]
        .notna()
        .sum()
    )

    if n_alt != 3074:
        raise RuntimeError(
            "Expected ALT available for "
            f"n=3074; got {n_alt}"
        )

    if n_log_alt != 3074:
        raise RuntimeError(
            "Expected log ALT available for "
            f"n=3074; got {n_log_alt}"
        )

    raw_skew = (
        locked["LBXSATSI"]
        .skew()
    )

    log_skew = (
        locked["log_ALT"]
        .skew()
    )

    if not np.isclose(
        raw_skew,
        14.901250,
        atol=1e-5,
    ):
        raise RuntimeError(
            "Raw ALT skewness does not reproduce "
            f"frozen QC: {raw_skew}"
        )

    if not np.isclose(
        log_skew,
        0.699786,
        atol=1e-5,
    ):
        raise RuntimeError(
            "Log ALT skewness does not reproduce "
            f"frozen QC: {log_skew}"
        )

    # --------------------------------------------------------
    # Common ALT domain
    #
    # BOTH models use exactly the same participants.
    # --------------------------------------------------------

    common_domain = (
        domain
        & df["log_ALT"].notna()
    )

    n_common = int(
        common_domain.sum()
    )

    if n_common != 3074:
        raise RuntimeError(
            "Expected ALT common domain "
            f"n=3074; got {n_common}"
        )

    # --------------------------------------------------------
    # Fit both models on identical n=3074 sample
    # --------------------------------------------------------

    rows = []

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary model, common sample",
            add_alt=False,
        )
    )

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary + log ALT",
            add_alt=True,
        )
    )

    results = pd.DataFrame(
        rows
    )

    # --------------------------------------------------------
    # Quantify change in focal sex-interaction estimate
    # --------------------------------------------------------

    reference = results.iloc[0]
    adjusted = results.iloc[1]

    delta_contrast = (
        adjusted["male_minus_female"]
        - reference["male_minus_female"]
    )

    abs_ref = abs(
        reference["male_minus_female"]
    )

    abs_adjusted = abs(
        adjusted["male_minus_female"]
    )

    percent_change_magnitude = (
        100
        * (abs_adjusted - abs_ref)
        / abs_ref
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
        out / "hepatic_sensitivity.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Display
    # --------------------------------------------------------

    print()
    print(
        "STAGE 2E: HEPATIC/METABOLIC BIOCHEMISTRY SENSITIVITY"
    )
    print("=" * 115)

    print(
        "Locked primary analytic domain: "
        f"{n_locked:,}"
    )

    print(
        "ALT common domain: "
        f"{n_common:,}"
    )

    print(
        "Excluded for missing ALT: "
        f"{n_locked - n_common:,}"
    )

    print()
    print(
        "Frozen Stage 2E specification: "
        "log-transformed ALT (log LBXSATSI)"
    )

    print(
        "Specification was frozen before "
        "outcome-model fitting."
    )

    print()
    print(
        "Both models below use exactly the same "
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
        "CHANGE AFTER LOG-ALT ADJUSTMENT"
    )
    print("-" * 115)

    print(
        "Male-minus-female contrast before "
        "ALT adjustment: "
        f"{reference['male_minus_female']:.6f}"
    )

    print(
        "Male-minus-female contrast after "
        "ALT adjustment:  "
        f"{adjusted['male_minus_female']:.6f}"
    )

    print(
        "Change in contrast:                 "
        f"{delta_contrast:.6f}"
    )

    print(
        "Percent change in absolute "
        "contrast magnitude: "
        f"{percent_change_magnitude:.2f}%"
    )

    print()
    print("SEX-SPECIFIC SLOPES")
    print("-" * 115)

    print(
        "Female before ALT adjustment: "
        f"{reference['female_beta']:.6f}"
    )

    print(
        "Female after ALT adjustment:  "
        f"{adjusted['female_beta']:.6f}"
    )

    print()

    print(
        "Male before ALT adjustment:   "
        f"{reference['male_beta']:.6f}"
    )

    print(
        "Male after ALT adjustment:    "
        f"{adjusted['male_beta']:.6f}"
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
        "ALT was modeled continuously on the "
        "natural-log scale."
    )

    print(
        "No participants were excluded on the basis "
        "of high or low ALT values."
    )

    print()


if __name__ == "__main__":
    main()
