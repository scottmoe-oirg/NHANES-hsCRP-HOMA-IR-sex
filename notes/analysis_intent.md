# NHANES Metabolic Follow-up: Initial Analysis Intent

Date: September 16, 2026

## Motivation

This project is an exploratory follow-up to the submitted hsCRP-HOMA-IR NHANES study.

During the methodological audit of Study #1, exploratory analyses suggested that sex modification of the HOMA-IR–hsCRP association might be important. Sex was not a primary effect modifier in Study #1, and this observation was not elevated to a primary finding.

## Initial Question

Among U.S. adults aged 20 years or older without diagnosed diabetes and with HbA1c <5.7%, does the adjusted association between HOMA-IR and systemic inflammation, measured by hsCRP, differ by sex?

## Initial Analysis

The first reconnaissance model will use the same final analytic population and survey-design principles as Study #1 (final analytic n = 3,079).

Primary reconnaissance specification:

    log(hsCRP) ~ HOMA_IR * sex + HbA1c + waist + age + nonHDL

The parameter of principal interest is the HOMA_IR × sex interaction.

Planned initial supporting analyses:

1. Survey-weighted descriptive characteristics by sex.
2. Sex-specific adjusted HOMA-IR slopes derived from the interaction model.
3. Raw hsCRP outcome as an outcome-scale sensitivity analysis.

A HOMA-IR × sex × age analysis will not be treated as part of the initial test. It will be considered subsequently if the two-way sex analysis provides sufficient scientific motivation.

## Project Boundary

Study #1 is frozen as submitted to medRxiv. Files in the original project will be treated as source/reference material and will not be modified by this follow-up project.
