# Stage 1 — Sex-Interaction Scale Robustness and Model Diagnostics

**Study:** NHANES metabolic follow-up  
**Analytic sample:** n = 3,079  
**Design denominator df:** 40  
**Primary scientific question:** Among U.S. adults age >=20 without diagnosed diabetes and with HbA1c <5.7%, does the adjusted association between HOMA-IR and hsCRP differ by sex?

## 1. Starting model

Initial reconnaissance model:

    log(hsCRP) ~ HOMA-IR * sex
                 + period * (HbA1c + waist + age + non-HDL)

Sex coding: female = 0, male = 1.

Initial result:

- HOMA-IR x sex: F(1,40) = 9.053325, p = 0.004521562
- Female HOMA-IR slope = 0.062864
  - SE = 0.018001
  - 95% CI = 0.026483 to 0.099245
  - p = 0.001183557
- Male HOMA-IR slope = 0.008172
  - SE = 0.011743
  - 95% CI = -0.015561 to 0.031905
  - p = 0.490491651
- Male minus female slope difference = -0.054692
  - 95% CI = -0.091428 to -0.017955

The formal interaction, rather than the separate within-sex p-values, is the relevant test of sex heterogeneity.

## 2. Pre-specified 2 x 2 scale sensitivity

All four combinations of raw/log outcome and raw/log exposure were fit with otherwise identical model structure.

| Outcome / exposure | Interaction F(1,40) | p | Female slope | Male slope |
|---|---:|---:|---:|---:|
| log hsCRP / raw HOMA-IR | 9.053325 | 0.004522 | 0.062864 | 0.008172 |
| raw hsCRP / raw HOMA-IR | 11.283858 | 0.001727 | 0.425004 | -0.033908 |
| log hsCRP / log HOMA-IR | 10.876567 | 0.002050 | 0.238919 | 0.016592 |
| raw hsCRP / log HOMA-IR | 15.311641 | 0.000345 | 1.282906 | -0.397836 |

**Stage 1 scale conclusion:** the direction and qualitative pattern of sex heterogeneity were preserved under all four parameterizations. The finding is therefore not dependent on one raw/log combination.

No model was selected on the basis of the smallest p-value.

## 3. Residual diagnostics: outcome scale

Overall residual summaries:

| Model | Residual SD | Skewness | Excess kurtosis | Max absolute residual |
|---|---:|---:|---:|---:|
| log hsCRP / raw HOMA-IR | 1.0420 | 0.3822 | 0.7937 | 4.5304 |
| raw hsCRP / raw HOMA-IR | 6.2884 | 11.8421 | 254.6930 | 180.4166 |
| log hsCRP / log HOMA-IR | 1.0418 | 0.3798 | 0.8060 | 4.5083 |
| raw hsCRP / log HOMA-IR | 6.2932 | 11.8145 | 253.9764 | 180.4489 |

The raw-hsCRP models showed extreme right-tail residual behavior. Logging hsCRP produced a major improvement in residual symmetry and tail behavior regardless of HOMA-IR scale.

Residual-vs-fitted plots supported the same conclusion: raw hsCRP produced a highly asymmetric residual geometry dominated by a small number of very large positive residuals, whereas log hsCRP produced a much more compact and approximately stable residual cloud.

**Outcome-scale conclusion:** log hsCRP is strongly favored as the primary outcome on model-specification grounds. This choice is not motivated by statistical significance; the raw-hsCRP interaction p-values were actually smaller.

## 4. Residual diagnostics: exposure scale

With log hsCRP fixed as the outcome, raw and log HOMA-IR produced nearly identical overall residual SD, skewness, and kurtosis.

Residual-vs-exposure plots showed that raw HOMA-IR is highly right-skewed, with most observations concentrated at low values and a sparse upper tail extending to approximately 39. Logging HOMA-IR spreads the dense lower range and compresses the upper tail.

However, neither raw nor log HOMA-IR produced a compelling gross residual pattern demonstrating obvious functional-form failure.

Binned residual means by pooled HOMA-IR quantile and sex fluctuated around zero under both exposure parameterizations. The plots did not provide strong evidence that one exposure scale was clearly superior.

The binned error bars were descriptive ordinary standard errors, not design-based inferential confidence intervals, and were used only as model-diagnostic aids.

## 5. HOMA-IR distribution and upper-tail influence

Observed HOMA-IR distribution in the analytic sample:

- 50th percentile: 1.991457
- 75th percentile: 3.218049
- 90th percentile: 4.962064
- 95th percentile: 6.662040
- 97.5th percentile: 8.437815
- 99th percentile: 11.311570
- maximum: 38.974321

Diagnostic upper-tail restrictions were then applied. These are sensitivity perturbations, not proposed participant exclusions.

### Raw HOMA-IR

| Restriction | n | Interaction p | Female slope | Male slope | Male - female |
|---|---:|---:|---:|---:|---:|
| Full sample | 3079 | 0.004522 | 0.062864 | 0.008172 | -0.054692 |
| Below 99th percentile | 3048 | 0.013336 | 0.095393 | 0.021941 | -0.073453 |
| Below 97.5th percentile | 3002 | 0.001625 | 0.118833 | 0.015282 | -0.103551 |
| Below 95th percentile | 2925 | 0.009464 | 0.121024 | 0.020917 | -0.100107 |

The interaction did not disappear when the sparse upper tail was removed. The estimated female slope and female-male separation became larger under the raw-HOMA restrictions. Thus the extreme upper tail is not generating the interaction; if anything, it attenuates the fitted raw-scale female slope.

### Log HOMA-IR

| Restriction | n | Interaction p | Female slope | Male slope | Male - female |
|---|---:|---:|---:|---:|---:|
| Full sample | 3079 | 0.002050 | 0.238919 | 0.016592 | -0.222327 |
| Below 99th percentile | 3048 | 0.004497 | 0.244022 | 0.018937 | -0.225085 |
| Below 97.5th percentile | 3002 | 0.001350 | 0.261160 | 0.000295 | -0.260865 |
| Below 95th percentile | 2925 | 0.004835 | 0.243248 | 0.000246 | -0.243002 |

The log-HOMA slopes and sex contrast were comparatively stable across upper-tail perturbations.

**Upper-tail conclusion:** the sex interaction is not an artifact of a few extremely high HOMA-IR observations. Log HOMA-IR shows greater coefficient stability to upper-tail restriction, but the residual diagnostics do not establish that log HOMA-IR is uniquely correct.

## 6. Survey-regression interpretation

Inference uses the NHANES survey design and Taylor-linearized covariance estimates. Therefore, textbook IID/homoscedastic residual assumptions are not being used as the basis of the ordinary OLS covariance formula.

Residual diagnostics remain important for a different reason: model specification. They inform us about functional form, extreme observations, variance structure, and whether the transformation yields a sensible conditional-mean model.

Design-based covariance and residual diagnostics therefore answer complementary questions:

- **Survey design:** how should uncertainty in the fitted coefficients be quantified?
- **Residual diagnostics:** does the proposed regression form adequately represent the observed relationship?

Good design-based standard errors do not rescue a poorly specified conditional-mean model, and well-behaved residuals do not justify ignoring the complex survey design.

## 7. Stage 1 conclusions frozen before Stage 2

1. Evidence of HOMA-IR x sex heterogeneity is present under all four raw/log outcome-exposure combinations tested.
2. Log hsCRP is strongly favored as the primary outcome because it dramatically improves residual behavior.
3. This outcome transformation is not being selected to obtain significance; the interaction is also present, and numerically more significant, on the raw hsCRP scale.
4. Raw versus log HOMA-IR is not decisively settled by residual diagnostics.
5. Log HOMA-IR gives notably more stable sex-specific slope estimates under progressive upper-tail restriction.
6. The sex interaction is not generated by the sparse extreme HOMA-IR tail.
7. Upper-tail restrictions are diagnostic sensitivities only and are not new analytic exclusion criteria.
8. These analyses establish robustness and model behavior; they do not establish mechanism, causality, or biological mediation.
9. No three-way HOMA-IR x sex x age interaction has yet been pursued.
10. The next planned scientific stage is to examine whether the observed sex heterogeneity depends materially on how adiposity is represented.

**Decision point:** proceed to Stage 2 without changing the Stage 1 interpretation in response to later findings.
