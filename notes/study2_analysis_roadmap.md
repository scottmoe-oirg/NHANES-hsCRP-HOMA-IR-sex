# NHANES Metabolic Follow-up: Study 2 Analysis Roadmap

**Date:** September 17, 2026

## 1. Study Question

Among U.S. adults aged ≥20 years without diagnosed diabetes and with
HbA1c <5.7%, does the adjusted association between HOMA-IR and systemic
inflammation (hsCRP) differ by sex?

This study follows the submitted hsCRP–HOMA-IR NHANES study but is a
separate analysis. Study 1 is frozen and will not be modified.

---

## 2. Analytic Population and Survey Design

Study 2 begins with the same final analytic population and survey
principles as Study 1:

- NHANES 2015–March 2020 pre-pandemic data
- Age ≥20 years
- No self-reported diagnosed diabetes
- HbA1c <5.7%
- Qualifying fasting-subsample weights
- Valid fasting glucose and insulin
- Positive HOMA-IR and hsCRP
- Complete primary-model covariates
- Final analytic n = 3,079
- NHANES fasting weights, strata, and PSU structure retained
- Taylor-linearized survey variance estimation
- Design denominator df = 40 in the initial analysis

The Study 1 pipeline has been independently reproduced in the Study 2
workspace and passed all locked validation targets before new analyses
were performed.

---

## 3. Primary Reconnaissance Model

The initial sex-modification model was specified before its result was
examined:

    log(hsCRP) ~ HOMA_IR * is_male
                 + C(period) * (HbA1c + waist + age + nonHDL)

Internal variable names:

    log_hs_CRP ~ HOMA_IR * is_male
                 + C(period) * (LBXGH + waist + AGE + Non_HDL)

Sex coding:

- Female = 0
- Male = 1

Therefore:

- HOMA_IR coefficient = adjusted female HOMA-IR slope
- HOMA_IR:is_male coefficient = male slope − female slope
- The formal HOMA-IR × sex interaction is the principal parameter of
  interest.

---

## 4. Initial Result

Analytic population:

- Female: n = 1,638; weighted population = 74,169,596 (51.39%)
- Male: n = 1,441; weighted population = 70,154,211 (48.61%)

Formal interaction:

    F(1,40) = 9.053325
    p = 0.004521562

Adjusted HOMA-IR slopes:

Female:

    beta = 0.062864
    SE = 0.018001
    95% CI = 0.026483 to 0.099245
    p = 0.001184

Male:

    beta = 0.008172
    SE = 0.011743
    95% CI = -0.015561 to 0.031905
    p = 0.490492

Estimated male–female difference in slopes:

    approximately -0.05469

This result provides initial evidence that the adjusted HOMA-IR–hsCRP
association differs by sex.

The interpretation is based on the formal interaction test, not on the
fact that one sex-specific slope is statistically significant and the
other is not.

---

## 5. Stage 1: Scale Robustness

Study 1 demonstrated that interaction results can be sensitive to the
scale on which HOMA-IR and hsCRP are modeled. Therefore, the first task
is to characterize the sex interaction under four reasonable
parameterizations.

Planned 2 × 2 analysis:

| Outcome | Exposure | Status |
|---|---|---|
| log(hsCRP) | raw HOMA-IR | Completed — primary reconnaissance |
| raw hsCRP | raw HOMA-IR | Planned |
| log(hsCRP) | log HOMA-IR | Planned |
| raw hsCRP | log HOMA-IR | Planned |

For each model record:

- analytic n
- design df
- formal HOMA-IR × sex interaction F statistic
- interaction p-value
- female slope, SE, 95% CI, and p-value
- male slope, SE, 95% CI, and p-value

These analyses are intended to assess scale dependence, not to select
the parameterization producing the smallest p-value.

No age-stratified or three-way interaction analysis will be performed
before this scale audit is examined.

---

## 6. Stage 2: Influence, Shape, and Adiposity

If Stage 1 supports continued investigation, characterize whether the
observed sex heterogeneity reflects a broad relationship or is
sensitive to influential observations, distributional tails, or model
shape.

### 6.1 Influence and Functional Form

Potential analyses include:

- visualization of HOMA-IR and hsCRP distributions by sex
- visualization of adjusted/predicted relationships
- assessment of influential observations
- restricted-range or tail sensitivity analyses where scientifically
  justified
- assessment of possible nonlinear HOMA-IR relationships

These analyses should characterize the structure of the data rather
than remove observations merely because they weaken an interaction.

### 6.2 Adiposity Representation

Evaluate how the estimated sex interaction changes under different
representations of adiposity.

Planned conceptual sequence:

    no adiposity adjustment
        ↓
    BMI
        ↓
    waist circumference
        ↓
    waist-to-height or other justified adiposity representation

The purpose is to ask whether differences in adiposity representation
account for, attenuate, strengthen, or otherwise alter the observed
sex difference in metabolic–inflammatory coupling.

No specification will be selected solely because it produces a more
favorable p-value.

---

## 7. Stage 3: Relationship Between Sex and Age Heterogeneity

Study 1 identified age heterogeneity in the HOMA-IR–hsCRP association.
Study 2 has now produced preliminary evidence of sex heterogeneity.

A later conditional question is:

> Are the age and sex patterns in metabolic–inflammatory coupling
> related?

This question will not initially be approached by automatically fitting
a HOMA-IR × sex × age interaction.

If justified by Stages 1 and 2, first characterize and visualize
sex-specific HOMA-IR–hsCRP relationships across adulthood.

Possible later analyses may include:

- sex-specific slopes within the previously defined age groups
- age modeled continuously where appropriate
- visualization of age-dependent sex differences
- formal higher-order interaction only if scientifically and
  statistically justified

The limited survey design degrees of freedom and the risk of unstable
subgroup estimates must be considered.

---

## 8. Conditional Physiological Investigation

If robust sex and/or age-by-sex structure remains after the preceding
analyses, investigate biologically plausible contributors.

Potential domains include:

- central versus general adiposity
- reproductive aging
- hormonal environment where NHANES measurements permit
- smoking
- statin or other medication use
- renal function
- liver/metabolic variables
- lipid physiology

These analyses are exploratory unless separately specified in advance.

Age must not be treated as equivalent to menopausal status, and
observed statistical patterns must not be interpreted as evidence for
a specific hormonal mechanism without appropriate data.

---

## 9. Interpretation Boundaries

This is a cross-sectional NHANES study.

Accordingly:

- association does not establish causation
- HOMA-IR is a surrogate measure of insulin resistance
- hsCRP is a systemic inflammatory marker, not a specific inflammatory
  pathway
- sex modification does not by itself identify a biological mechanism
- statistical interaction may depend on variable scale
- subgroup-specific significance does not establish between-group
  heterogeneity; the interaction must be tested directly
- sensitivity analyses that weaken the initial result are scientifically
  informative and will be reported rather than treated as failed
  analyses

---

## 10. Broader Scientific Motivation

The broader program asks how much physiological heterogeneity exists
within populations classified as clinically normal.

HbA1c <5.7% defines a population without conventional laboratory
evidence of prediabetes, but does not imply metabolic homogeneity.

Study 1 investigated hidden structure by age.

Study 2 investigates whether metabolic–inflammatory coupling also
differs by sex.

The broader conceptual theme is:

> Study physiology before diagnosis.

---

## 11. Current Decision Point

As of September 17, 2026, the next analysis is fixed:

**Complete the remaining three cells of the HOMA-IR/hsCRP scale
sensitivity matrix before pursuing age stratification, higher-order
interactions, or explanatory physiological analyses.**

The result of that scale audit will determine the next stage of the
study.
