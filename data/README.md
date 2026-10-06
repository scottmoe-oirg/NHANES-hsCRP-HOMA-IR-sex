# NHANES Raw Data

This directory contains the raw public-use NHANES data files required to
reproduce the analyses for the study:

**Sex Differences in Metabolic–Inflammatory Coupling Below the HbA1c
Prediabetes Threshold: NHANES 2015–March 2020**

The raw NHANES `.xpt` files are not distributed with this repository.
They are publicly available from the National Center for Health Statistics
(NCHS) and should be downloaded from the official NHANES website.

Place all downloaded files in:

    data/raw/

The analysis scripts reconstruct the analytic dataset directly from these
public-use files.

## NHANES periods

The study combines:

- NHANES 2015–2016
- NHANES 2017–March 2020 pre-pandemic

NHANES filenames ending in `_I` correspond to the 2015–2016 cycle.
Files beginning with `P_` correspond to the 2017–March 2020
pre-pandemic release.

## Core files required for the validated analytic cohort

The following 20 files are read directly by `src/data.py` and are required
to reconstruct the locked analytic population and the variables used by the
core analysis framework.

### NHANES 2015–2016

    DEMO_I.xpt
    BMX_I.xpt
    DIQ_I.xpt
    FASTQX_I.xpt
    GHB_I.xpt
    GLU_I.xpt
    INS_I.xpt
    HSCRP_I.xpt
    TCHOL_I.xpt
    HDL_I.xpt

### NHANES 2017–March 2020 pre-pandemic

    P_DEMO.xpt
    P_BMX.xpt
    P_DIQ.xpt
    P_FASTQX.xpt
    P_GHB.xpt
    P_GLU.xpt
    P_INS.xpt
    P_HSCRP.xpt
    P_TCHOL.xpt
    P_HDL.xpt

## Additional files used for sensitivity analyses and diagnostics

The following files support the targeted smoking, alcohol, medication,
renal, hepatic, pregnancy/reproductive-health, and complete-blood-count
analyses included in the archived workflow.

### NHANES 2015–2016

    ALQ_I.xpt
    BIOPRO_I.xpt
    CBC_I.xpt
    RHQ_I.xpt
    RXQ_RX_I.xpt
    SMQ_I.xpt

### NHANES 2017–March 2020 pre-pandemic

    P_ALQ.xpt
    P_BIOPRO.xpt
    P_CBC.xpt
    P_RHQ.xpt
    P_RXQ_RX.xpt
    P_SMQ.xpt

### Prescription-drug reference file

    RXQ_DRUG.xpt


## Role of the supplemental files

The primary HOMA-IR × sex model is constructed from the demographic,
examination, fasting, glycemic, lipid, and inflammatory data required by
the locked analytic framework.

Additional NHANES components support the targeted sensitivity analyses and
descriptive diagnostics reported with the study. These include:

- alcohol use
- smoking
- prescription medication and statin use
- renal function
- hepatic biomarkers
- pregnancy status
- reproductive-health quality control
- complete blood count and leukocyte differential measures

Not every file listed above is therefore required to fit the primary model
alone. The full collection is required to reproduce the complete set of analyses, 
sensitivity analyses, and diagnostics retained in the public study repository.” 

## Data privacy

Only publicly available, de-identified NHANES public-use data are used.
No restricted-use NHANES data are required.

Raw NHANES files and participant-level derived datasets are intentionally
excluded from version control. Aggregate analysis outputs needed to verify
the reported results are retained in the repository.
