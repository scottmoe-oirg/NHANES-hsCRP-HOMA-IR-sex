#!/usr/bin/env python3

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


STATIN_CLASS = "HMG-COA REDUCTASE INHIBITORS (STATINS)"


# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

def load_config(path: Path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def clean_text(x):
    """
    Convert SAS bytes/strings to clean uppercase text.
    """

    if pd.isna(x):
        return ""

    if isinstance(x, bytes):
        x = x.decode(
            "utf-8",
            errors="ignore",
        )

    return str(x).strip().upper()


def build_statin_status(
    medication_dir: Path,
):
    """
    Construct the pre-specified participant-level statin exposure.

    Definite statin:
        At least one reported prescription medication contains
        an ingredient whose NHANES ingredient-level therapeutic
        class is HMG-COA REDUCTASE INHIBITORS (STATINS).

    Definite nonstatin:
        No prescription medication use, OR prescription medication
        use with all reported medications classifiable and no
        statin detected.

    Ambiguous:
        Reported prescription medication use with >=1 unidentified
        medication and no detected statin, OR refused/don't know
        prescription-use response.

    Participants with an identified statin remain definite statin
    users even if another reported medication is unidentified.
    """

    drug_path = (
        medication_dir
        / "RXQ_DRUG.xpt"
    )

    rx_paths = [
        medication_dir / "RXQ_RX_I.xpt",
        medication_dir / "P_RXQ_RX.xpt",
    ]

    if not drug_path.exists():
        raise FileNotFoundError(
            f"Missing drug classification file: "
            f"{drug_path}"
        )

    for path in rx_paths:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing prescription file: "
                f"{path}"
            )

    # --------------------------------------------------------
    # Drug classification file
    # --------------------------------------------------------

    drug = pd.read_sas(
        drug_path,
        format="xport",
    )

    drug["RXDDRGID_clean"] = (
        drug["RXDDRGID"]
        .map(clean_text)
    )

    ingredient_class_cols = [
        col
        for col in drug.columns
        if col.startswith("RXDICN")
    ]

    if len(ingredient_class_cols) != 18:
        raise RuntimeError(
            "Expected 18 ingredient-level therapeutic "
            f"class columns; got {len(ingredient_class_cols)}"
        )

    def contains_statin(row):
        return any(
            clean_text(row[col])
            == STATIN_CLASS
            for col in ingredient_class_cols
        )

    drug["contains_statin"] = (
        drug.apply(
            contains_statin,
            axis=1,
        )
    )

    n_statin_ids = int(
        drug["contains_statin"].sum()
    )

    if n_statin_ids != 14:
        raise RuntimeError(
            "Expected 14 statin-containing drug IDs "
            f"from frozen QC; got {n_statin_ids}"
        )

    drug_map = (
        drug[
            [
                "RXDDRGID_clean",
                "contains_statin",
            ]
        ]
        .set_index("RXDDRGID_clean")[
            "contains_statin"
        ]
        .to_dict()
    )

    # --------------------------------------------------------
    # Prescription medication files
    # --------------------------------------------------------

    frames = []

    for path in rx_paths:

        rx_part = pd.read_sas(
            path,
            format="xport",
        )

        needed = [
            "SEQN",
            "RXDUSE",
            "RXDDRUG",
            "RXDDRGID",
        ]

        missing = [
            x
            for x in needed
            if x not in rx_part.columns
        ]

        if missing:
            raise RuntimeError(
                f"{path.name} missing columns: "
                f"{missing}"
            )

        frames.append(
            rx_part[needed].copy()
        )

    rx = pd.concat(
        frames,
        ignore_index=True,
    )

    rx["RXDDRGID_clean"] = (
        rx["RXDDRGID"]
        .map(clean_text)
    )

    # --------------------------------------------------------
    # Frozen linkage QC
    # --------------------------------------------------------

    meaningful = (
        rx["RXDDRGID_clean"] != ""
    )

    n_nonblank = int(
        meaningful.sum()
    )

    n_matched = int(
        rx.loc[
            meaningful,
            "RXDDRGID_clean",
        ]
        .isin(drug_map)
        .sum()
    )

    if n_nonblank != 37612:
        raise RuntimeError(
            "Expected 37,612 prescription rows with "
            f"nonblank RXDDRGID; got {n_nonblank}"
        )

    if n_matched != 37612:
        raise RuntimeError(
            "Expected all 37,612 nonblank drug IDs "
            f"to map; got {n_matched}"
        )

    # --------------------------------------------------------
    # Row-level classification
    # --------------------------------------------------------

    rx["mapped_drug"] = (
        rx["RXDDRGID_clean"]
        .isin(drug_map)
    )

    rx["statin_row"] = (
        rx["RXDDRGID_clean"]
        .map(drug_map)
        .fillna(False)
        .astype(bool)
    )

    rx["unidentified_med_row"] = (
        (rx["RXDUSE"] == 1)
        & (~rx["mapped_drug"])
    )

    # --------------------------------------------------------
    # Collapse to participant level
    # --------------------------------------------------------

    person_rows = []

    for seqn, g in rx.groupby("SEQN"):

        rxduse_values = set(
            g["RXDUSE"]
            .dropna()
            .tolist()
        )

        any_rx = (
            1.0 in rxduse_values
        )

        no_rx = (
            rxduse_values == {2.0}
        )

        refused_unknown = any(
            value in {7.0, 9.0}
            for value in rxduse_values
        )

        statin = bool(
            g["statin_row"].any()
        )

        unidentified = bool(
            g["unidentified_med_row"].any()
        )

        if statin:
            status = "definite_statin"

        elif (
            no_rx
            and not refused_unknown
        ):
            status = "definite_nonstatin"

        elif (
            any_rx
            and not unidentified
            and not refused_unknown
        ):
            status = "definite_nonstatin"

        else:
            status = "ambiguous"

        person_rows.append(
            {
                "SEQN": seqn,
                "statin_status": status,
            }
        )

    person = pd.DataFrame(
        person_rows
    )

    person["statin_use"] = np.where(
        person["statin_status"]
        == "definite_statin",
        1.0,
        np.where(
            person["statin_status"]
            == "definite_nonstatin",
            0.0,
            np.nan,
        ),
    )

    return person


def fit_model(
    df,
    domain,
    label,
    add_statin=False,
):

    base = (
        "log_hs_CRP ~ HOMA_IR * is_male "
        "+ C(period) * "
        "(LBXGH + waist + AGE + Non_HDL)"
    )

    if add_statin:
        formula = (
            base
            + " + statin_use"
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
        "--medication-dir",
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
            / "sex_statin_sensitivity"
        ),
    )

    args = parser.parse_args()

    config = load_config(
        Path(args.config)
    )

    # --------------------------------------------------------
    # Build validated Study 2 dataset
    # --------------------------------------------------------

    df, domain, flow = (
        build_combined_dataset(
            Path(args.data_dir),
            config,
        )
    )

    if int(domain.sum()) != 3079:
        raise RuntimeError(
            "Expected primary analytic domain "
            f"n=3079; got {int(domain.sum())}"
        )

    # --------------------------------------------------------
    # Construct frozen statin exposure
    # --------------------------------------------------------

    statin = build_statin_status(
        Path(args.medication_dir)
    )

    # --------------------------------------------------------
    # Merge onto validated dataset
    # --------------------------------------------------------

    original_rows = len(df)

    df = df.merge(
        statin[
            [
                "SEQN",
                "statin_status",
                "statin_use",
            ]
        ],
        on="SEQN",
        how="left",
        validate="one_to_one",
        sort=False,
    )

    if len(df) != original_rows:
        raise RuntimeError(
            "Statin merge changed dataset row count."
        )

    # --------------------------------------------------------
    # Frozen QC counts within locked domain
    # --------------------------------------------------------

    locked = df.loc[
        domain
    ].copy()

    counts = (
        locked["statin_status"]
        .value_counts(
            dropna=False
        )
    )

    expected_counts = {
        "definite_nonstatin": 2755,
        "definite_statin": 286,
        "ambiguous": 38,
    }

    for status, expected in expected_counts.items():

        observed = int(
            counts.get(
                status,
                0,
            )
        )

        if observed != expected:
            raise RuntimeError(
                f"Expected {status} n={expected}; "
                f"got {observed}"
            )

    # --------------------------------------------------------
    # Common classifiable domain
    #
    # Both models use exactly the same participants.
    # --------------------------------------------------------

    statin_complete = (
        df["statin_use"].notna()
    )

    common_domain = (
        domain
        & statin_complete
    )

    n_common = int(
        common_domain.sum()
    )

    if n_common != 3041:
        raise RuntimeError(
            "Expected statin-classifiable domain "
            f"n=3041; got {n_common}"
        )

    # --------------------------------------------------------
    # Fit BOTH models on identical n=3041 sample
    # --------------------------------------------------------

    rows = []

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary model, common sample",
            add_statin=False,
        )
    )

    rows.append(
        fit_model(
            df,
            common_domain,
            "Primary + statin use",
            add_statin=True,
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
        out / "statin_sensitivity.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Display
    # --------------------------------------------------------

    print()
    print(
        "STAGE 2C: STATIN ADJUSTMENT SENSITIVITY"
    )
    print("=" * 110)

    print(
        "Locked primary analytic domain: "
        f"{int(domain.sum()):,}"
    )

    print(
        "Statin-classifiable common domain: "
        f"{n_common:,}"
    )

    print(
        "Excluded for ambiguous statin status: "
        f"{int(domain.sum()) - n_common:,}"
    )

    print()
    print(
        "Frozen classification within "
        "locked primary domain:"
    )

    print(
        "  definite statin:    286"
    )

    print(
        "  definite nonstatin: 2,755"
    )

    print(
        "  ambiguous:           38"
    )

    print()
    print(
        "Both models below use exactly the same "
        f"n={n_common:,} participants."
    )

    print()
    print("RESULTS")
    print("-" * 110)

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
        "CHANGE AFTER STATIN ADJUSTMENT"
    )
    print("-" * 110)

    print(
        "Male-minus-female contrast before "
        "statin adjustment: "
        f"{reference['male_minus_female']:.6f}"
    )

    print(
        "Male-minus-female contrast after "
        "statin adjustment:  "
        f"{adjusted['male_minus_female']:.6f}"
    )

    print(
        "Absolute change in contrast:        "
        f"{delta_contrast:.6f}"
    )

    print(
        "Percent change in absolute "
        "contrast magnitude: "
        f"{percent_change_magnitude:.2f}%"
    )

    print()
    print("SEX-SPECIFIC SLOPES")
    print("-" * 110)

    print("Female:")

    print(
        "  before statin adjustment = "
        f"{reference['female_beta']:.6f}"
    )

    print(
        "  after statin adjustment  = "
        f"{adjusted['female_beta']:.6f}"
    )

    print()

    print("Male:")

    print(
        "  before statin adjustment = "
        f"{reference['male_beta']:.6f}"
    )

    print(
        "  after statin adjustment  = "
        f"{adjusted['male_beta']:.6f}"
    )

    print()
    print("MODELS")
    print("-" * 110)

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
        "Statin definition and handling of "
        "ambiguous medication records were "
        "frozen before fitting these models."
    )

    print(
        "Statin use = at least one reported "
        "prescription medication containing an "
        "ingredient classified by NHANES as "
        "HMG-COA reductase inhibitors (statins)."
    )

    print(
        "Participants with unresolved statin "
        "status were excluded from both models."
    )

    print()


if __name__ == "__main__":
    main()
