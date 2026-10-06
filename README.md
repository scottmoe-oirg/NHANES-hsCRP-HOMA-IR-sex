# Sex Differences in Metabolic–Inflammatory Coupling Below the HbA1c Prediabetes Threshold

This repository contains the computational analysis, reproducibility materials,
manuscript, and supplementary material for the study:

**Sex Differences in Metabolic–Inflammatory Coupling Below the HbA1c
Prediabetes Threshold: NHANES 2015–March 2020**

The study examines whether the adjusted association between insulin resistance,
estimated using the homeostatic model assessment of insulin resistance
(HOMA-IR), and systemic inflammation, measured using high-sensitivity
C-reactive protein (hsCRP), differs by sex among U.S. adults below the
HbA1c threshold for prediabetes.

The analysis uses publicly available data from the National Health and
Nutrition Examination Survey (NHANES), combining the 2015–2016 cycle with
the 2017–March 2020 pre-pandemic release.

## Primary analysis

The locked primary interaction model is

```text
log(hsCRP) ~ HOMA-IR * sex
             + period * (HbA1c + waist circumference + age + non-HDL cholesterol)
```

The analytic cohort includes adults aged 20 years or older without diagnosed
diabetes, with HbA1c <5.7%, qualifying fasting-subsample weights, valid fasting
glucose and insulin measurements, positive HOMA-IR and hsCRP values, and
complete primary-model covariates.

The final analytic sample contains **3,079 participants** and represents
approximately **144.3 million U.S. adults** under the combined NHANES survey
weights.

In the primary model, the HOMA-IR-by-sex interaction was:

```text
F(1,40) = 9.053325
p = 0.004522
```

Estimated HOMA-IR slopes on the log-hsCRP scale were:

```text
Women:  0.062864
Men:    0.008172
Male-minus-female contrast: -0.054692
95% CI: [-0.091428, -0.017955]
```

These estimates describe adjusted associations and should not be interpreted
as causal effects.

## Survey design

Analyses account for the NHANES complex survey design, including fasting
subsample weights, masked variance strata, and primary sampling units.

The 2015–2016 and 2017–March 2020 fasting weights are rescaled according to
their respective 2.0-year and 3.2-year contributions to the combined
5.2-year study period. Stratum and PSU identifiers are made unique across
periods.

Analytic-domain estimation retains the full qualifying survey-design frame
and sets score contributions to zero outside the analytic domain rather than
treating the analytic cohort as an independently sampled dataset.

## Repository structure

```text
config/       Analysis configuration
data/         Raw-data instructions; NHANES XPT files are not distributed
figures/      Study figures
manuscript/   Manuscript, supplement, bibliography, and build script
notes/        Analysis intent, decision history, and validation records
output/       Study #2 analysis and sensitivity-analysis outputs
outputs/      Core inherited-framework validation outputs
scripts/      Study #2 analysis, sensitivity, QC, and figure scripts
src/          Core data-construction and survey-regression code
validation/   Independent R validation of the primary sex-interaction model
```

## Raw data

Raw NHANES XPT files are not included in this repository. They are publicly
available from the National Center for Health Statistics.

Place the required files in:

```text
data/raw/
```

The exact 33-file dependency set, including the 20 files required for the
core analytic cohort and 13 files used for sensitivity analyses and
diagnostics, is documented in:

```text
data/README.md
```

No restricted-use NHANES data are required.

## Python environment

The analysis was conducted with:

```text
numpy==2.3.4
pandas==3.0.5
patsy==1.0.2
scipy==1.18.0
statsmodels==0.14.6
```

The complete Python dependency specification is provided in
`requirements.txt`.

A local environment can be created with, for example:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Reproducing the inherited analysis framework

The survey-analysis framework used for this study was inherited from a
preceding analysis of age-dependent HOMA-IR–hsCRP associations in the same
NHANES cohort.

To verify that the inherited computational framework reproduces the locked
preceding-study results:

```bash
python3 run_analysis.py --mode primary
```

The expected validation targets are stored in
`config/analysis_config.toml`, and the recorded reproduction check is in:

```text
notes/study1_reproduction_validation.json
```

This validation step is distinct from the primary sex-interaction analysis
reported in the present study.

## Reproducing the Study #2 primary analysis

Run:

```bash
python3 scripts/sex_reconnaissance.py
```

The expected primary result is:

```text
n = 3079
design df = 40
HOMA-IR × sex: F(1,40) = 9.053325, p = 0.004521562
```

Primary model outputs are written to:

```text
output/sex_reconnaissance/
```

## Sensitivity analyses and diagnostics

The `scripts/` directory contains the analysis programs used for the
study's targeted robustness and diagnostic analyses, including:

- exposure-scale and upper-tail sensitivity
- alternative adiposity specifications
- smoking adjustment
- statin-use adjustment
- renal-function adjustment
- hepatic-biomarker adjustment
- pregnancy exclusion
- alcohol adjustment
- restriction to hsCRP <=10 mg/L
- descriptive complete-blood-count assessment of the upper hsCRP group
- residual diagnostics and study figures

Aggregate outputs are retained under `output/`. Participant-level derived
datasets are excluded from version control.

## Independent R validation

The primary HOMA-IR-by-sex interaction model was independently reproduced
using the R `survey` package.

The validation script is:

```text
validation/validate_primary_sex_R.R
```

The independent implementation reproduced the locked primary interaction
test and sex-specific slopes to numerical tolerance.

This cross-software validation applies to the primary model; the complete
set of sensitivity analyses was not independently duplicated in R.

## Manuscript

The manuscript and supplementary material are provided in `manuscript/`.

To rebuild both PDFs:

```bash
./manuscript/compile_study2.sh
```

The build requires a LaTeX installation with BibTeX.

Canonical outputs are:

```text
manuscript/NHANES_sex_interaction_manuscript.pdf
manuscript/NHANES_sex_interaction_supplement.pdf
```

## Data privacy

This study uses only publicly available, de-identified NHANES public-use
data. Raw NHANES XPT files and participant-level derived datasets are not
distributed in this repository.

## Authors

**Scott Moerschbacher**  
Open Inquiry Research Group

**Janice Baker**  
Open Inquiry Research Group  
Department of Sociology, University of Nevada, Reno

## Citation

Citation information and the permanent Zenodo DOI will be added when the
public repository release is archived.

## License

License information will be added before the public release.
