Study #2 Analysis Log

## Study #2 Primary Interaction Model

### Scientific question

Study #2 evaluates whether the adjusted association between insulin resistance,
measured by HOMA-IR, and systemic inflammation, measured by hsCRP, differs by
sex among U.S. adults aged 20 years or older without diagnosed diabetes and
with HbA1c <5.7%.

The analysis uses the same NHANES 2015–March 2020 pre-pandemic fasting
subsample framework developed and validated in Study #1. The locked Study #2
analytic domain contains n=3,079 participants, with survey design denominator
degrees of freedom equal to 40.

The focal scientific question is one of effect modification: whether the
HOMA-IR–hsCRP slope differs between females and males. The relevant statistical
test is therefore the formal HOMA-IR × sex interaction rather than a comparison
of whether the sex-specific slopes are individually statistically significant.

### Locked primary model

The Study #2 primary interaction model is:

log(hsCRP) ~ HOMA_IR * is_male
             + C(period) * (LBXGH + waist + AGE + Non_HDL)

where:

- hsCRP is modeled on the natural-log scale.
- HOMA-IR is retained on its raw scale.
- is_male is coded 0 for female and 1 for male.
- LBXGH is HbA1c.
- waist is waist circumference.
- AGE is continuous age.
- Non_HDL is non-HDL cholesterol.
- period distinguishes the 2015–2016 and 2017–March 2020 NHANES periods.

The period-specific adjustment structure is retained from Study #1 so that the
associations of HbA1c, waist circumference, age, and non-HDL cholesterol with
the outcome are not forced to be identical across the two NHANES periods.

Because sex is the focal effect modifier in Study #2, is_male is not also
included as an ordinary adjustment covariate outside the HOMA-IR × sex
interaction.

The validated NHANES complex-survey analysis machinery from Study #1 is used
without modification, including fasting-subsample weights, strata, primary
sampling units, domain score-zeroing, Taylor-linearized covariance estimation,
and survey-design denominator degrees of freedom.

### Primary interaction result

In the locked n=3,079 analytic sample, the HOMA-IR × sex interaction was:

F(1,40) = 9.053325, p = 0.004522.

The adjusted HOMA-IR slope among females was:

beta = 0.062864
SE = 0.018001
95% CI: 0.026483 to 0.099245
p = 0.001184.

The adjusted HOMA-IR slope among males was:

beta = 0.008172
SE = 0.011743
95% CI: -0.015561 to 0.031905
p = 0.490492.

The male-minus-female HOMA-IR slope contrast was approximately:

-0.054692.

Because female is the reference category, the HOMA_IR coefficient represents
the female HOMA-IR slope and the HOMA_IR:is_male coefficient represents the
male-minus-female difference in slopes.

### Interpretation

The formal interaction provides evidence that the adjusted association between
HOMA-IR and log-hsCRP differs by sex under the Study #2 primary interaction
model. The estimated HOMA-IR–hsCRP association is more positive among females
than among males.

This result should not be described simply as "significant in females but not
significant in males." A difference in subgroup-specific statistical
significance does not itself demonstrate that the subgroup slopes differ. The
formal HOMA-IR × sex interaction directly tests that difference.

Because the outcome is log-transformed, the sex-specific slopes are additive
on the log-hsCRP scale. For descriptive interpretation on the original hsCRP
scale, exp(beta)-1 gives the proportional change in hsCRP associated with a
one-unit increase in HOMA-IR under the fitted model. The primary inferential
scale, however, remains log-hsCRP.

The primary result is an adjusted association and should not be interpreted as
evidence that insulin resistance causes inflammation, that inflammation causes
insulin resistance, or that sex itself identifies the biological mechanism
responsible for the difference in slopes.

Status: Study #2 primary interaction model locked.


## Stage 1 — Scale and Influence Sensitivity

### Rationale

The Study #2 primary interaction model uses raw HOMA-IR as the focal exposure
and log-transformed hsCRP as the outcome. Because HOMA-IR is right-skewed, the
primary interaction could potentially be sensitive to exposure scale or be
driven disproportionately by participants in the extreme upper tail of the
HOMA-IR distribution.

Stage 1 was therefore designed to evaluate whether the HOMA-IR × sex
interaction is stable to:

- log transformation of HOMA-IR;
- removal of progressively more of the extreme upper HOMA-IR tail; and
- examination of model residuals under the raw and log HOMA-IR
  specifications.

These analyses were used to evaluate model stability and influence rather than
to search for the specification producing the smallest p-value.

### HOMA-IR distribution

Within the locked Study #2 analytic domain, the HOMA-IR distribution was:

50th percentile: 1.991457  
75th percentile: 3.218049  
90th percentile: 4.962064  
95th percentile: 6.662040  
97.5th percentile: 8.437815  
99th percentile: 11.311570  
Maximum: 38.974321

The substantial difference between the upper percentiles and maximum confirmed
the presence of a long right tail and motivated explicit influence
sensitivities.

### Raw versus log HOMA-IR

In the full n=3,079 analytic sample, the primary raw-HOMA-IR model produced:

F(1,40) = 9.053325, p = 0.004522.

The female HOMA-IR slope was 0.062864, the male slope was 0.008172, and the
male-minus-female contrast was -0.054692.

When HOMA-IR was natural-log transformed, the interaction remained:

F(1,40) = 10.876567, p = 0.002050.

On the log-HOMA-IR scale, the female slope was 0.238919 and the male slope was
0.016592, giving a male-minus-female contrast of -0.222327.

The numerical slopes from the raw and log models are not directly comparable
in magnitude because they correspond to different exposure scales. The
important observation is that the estimated sex difference remained in the
same direction and the formal interaction persisted under both
parameterizations.

### Upper-tail influence sensitivity

To determine whether the interaction was being driven by a small number of
participants with extreme HOMA-IR values, the models were repeated after
excluding observations above the 99th, 97.5th, and 95th percentiles of the
HOMA-IR distribution.

#### Below the 99th percentile

Raw HOMA-IR, n=3,048:

F(1,40) = 6.706342, p = 0.013336.

Female slope = 0.095393.  
Male slope = 0.021941.  
Male-minus-female contrast = -0.073453.

Log HOMA-IR:

F(1,40) = 9.065661, p = 0.004497.

Female slope = 0.244022.  
Male slope = 0.018937.  
Male-minus-female contrast = -0.225085.

#### Below the 97.5th percentile

Raw HOMA-IR, n=3,002:

F(1,40) = 11.429335, p = 0.001625.

Female slope = 0.118833.  
Male slope = 0.015282.  
Male-minus-female contrast = -0.103551.

Log HOMA-IR:

F(1,40) = 11.875342, p = 0.001350.

Female slope = 0.261160.  
Male slope = 0.000295.  
Male-minus-female contrast = -0.260865.

#### Below the 95th percentile

Raw HOMA-IR, n=2,925:

F(1,40) = 7.431888, p = 0.009464.

Female slope = 0.121024.  
Male slope = 0.020917.  
Male-minus-female contrast = -0.100107.

Log HOMA-IR:

F(1,40) = 8.903035, p = 0.004835.

Female slope = 0.243248.  
Male slope = 0.000246.  
Male-minus-female contrast = -0.243002.

### Interpretation of the trimming analyses

The HOMA-IR × sex interaction persisted after removing the upper 1%, 2.5%,
and 5% of the HOMA-IR distribution and persisted under both raw and
log-transformed HOMA-IR specifications.

The interaction therefore does not appear to be generated solely by a small
number of participants with extremely high HOMA-IR values.

The fact that the raw-HOMA-IR sex contrast becomes larger in magnitude after
upper-tail trimming also argues against the interpretation that the observed
primary interaction exists merely because extreme HOMA-IR observations are
artificially producing the sex difference.

These trimmed analyses are sensitivity analyses rather than alternative
primary populations. The full locked analytic sample remains the primary
sample.

### Residual diagnostics

Residual diagnostics were examined for both the raw-HOMA-IR and
log-HOMA-IR specifications.

The principal features considered were whether residuals showed:

- systematic curvature or other structure across fitted values;
- strong changes in residual spread across the fitted range;
- a small number of observations dominating the residual pattern;
- substantial departures from the behavior expected under the fitted model;
  or
- a clear qualitative improvement after transforming HOMA-IR.

Both specifications produced acceptable residual behavior. Log transformation
of HOMA-IR produced at most a marginal improvement in some residual features;
it did not reveal that the raw-HOMA-IR specification was fundamentally
misspecified.

Residual diagnostics therefore did not provide a compelling reason to replace
the locked raw-HOMA-IR specification.

### Decision on HOMA-IR scale

Raw HOMA-IR was retained for the Study #2 primary interaction model.

This decision was based on several considerations:

1. The raw-HOMA-IR model produced acceptable residual diagnostics.
2. The HOMA-IR × sex interaction persisted after log transformation.
3. The interaction persisted after progressively removing the extreme upper
   HOMA-IR tail.
4. Raw HOMA-IR provides a straightforward interpretation in terms of a
   one-unit increase in HOMA-IR.
5. Retaining raw HOMA-IR preserves continuity with Study #1.
6. The modest differences between the raw and log specifications did not
   justify changing the primary model after inspection of the results.

The log-HOMA-IR analysis is retained as a transparent sensitivity analysis
rather than replacing the prespecified primary exposure scale.

### Overall Stage 1 interpretation

Stage 1 provides evidence that the estimated sex difference in the adjusted
HOMA-IR–hsCRP association is not simply an artifact of the raw HOMA-IR scale
or of a few observations in the extreme upper tail of the HOMA-IR
distribution.

Both raw and log HOMA-IR models support the same qualitative conclusion, and
progressively stronger upper-tail trimming does not eliminate the interaction.
Residual diagnostics do not identify a sufficiently important deficiency in
the raw-HOMA-IR specification to justify changing the primary model.

Accordingly, the Study #2 primary interaction model remains:

log(hsCRP) ~ HOMA_IR * is_male
             + C(period) * (LBXGH + waist + AGE + Non_HDL)

with raw HOMA-IR and log-transformed hsCRP.

Status: Stage 1 complete and frozen.

Stage 2 — Expanded Covariate / Adjustment Sensitivities

Stage 2A — Adiposity Sensitivity

Because adiposity is closely related to both insulin resistance and systemic inflammation, the sensitivity of the HOMA-IR × sex interaction to alternative representations of adiposity was examined. The locked Study #2 primary interaction model uses waist circumference as the primary adiposity adjustment. Alternative specifications using BMI and waist-to-height ratio were evaluated, together with a model containing no adiposity adjustment.

The four specifications were:

* No adiposity adjustment
* BMI adjustment
* Waist circumference adjustment (locked primary specification)
* Waist-to-height ratio adjustment

Both available-case and common-complete-case analyses were performed. BMI was available for 3,072 of the 3,079 participants in the locked analytic domain, waist circumference for all 3,079, height for 3,074, and waist-to-height ratio for 3,074. The common complete-case sample across all adiposity specifications contained n=3,072 participants.

Available-case results:

* No adiposity adjustment: n=3,079; F(1,40)=8.5058, p=0.00578; female HOMA-IR slope=0.18922; male slope=0.09376; male-minus-female contrast=-0.09546 (95% CI -0.16161 to -0.02931).
* BMI adjustment: n=3,072; F(1,40)=4.8638, p=0.03324; female slope=0.06061; male slope=0.01828; contrast=-0.04233 (95% CI -0.08113 to -0.00354).
* Waist circumference adjustment (primary): n=3,079; F(1,40)=9.0533, p=0.00452; female slope=0.06286; male slope=0.00817; contrast=-0.05469 (95% CI -0.09143 to -0.01796).
* Waist-to-height ratio adjustment: n=3,074; F(1,40)=5.0256, p=0.03059; female slope=0.05216; male slope=0.01057; contrast=-0.04159 (95% CI -0.07909 to -0.00410).

Common complete-case results (n=3,072 for all specifications):

* No adiposity adjustment: F(1,40)=8.5096, p=0.00577; female slope=0.18908; male slope=0.09360; contrast=-0.09548 (95% CI -0.16163 to -0.02933).
* BMI adjustment: F(1,40)=4.8638, p=0.03324; female slope=0.06061; male slope=0.01828; contrast=-0.04233 (95% CI -0.08113 to -0.00354).
* Waist circumference adjustment: F(1,40)=9.0537, p=0.00452; female slope=0.06297; male slope=0.00822; contrast=-0.05475 (95% CI -0.09152 to -0.01798).
* Waist-to-height ratio adjustment: F(1,40)=5.0363, p=0.03042; female slope=0.05221; male slope=0.01055; contrast=-0.04165 (95% CI -0.07917 to -0.00414).

Adiposity adjustment substantially attenuated the sex-specific HOMA-IR–hsCRP slopes relative to the model without adiposity adjustment. However, the estimated male-minus-female HOMA-IR slope contrast remained negative under BMI, waist circumference, and waist-to-height ratio adjustment, indicating that the adjusted HOMA-IR–hsCRP association remained more positive among females than males across these alternative adiposity representations.

Available-case and common-complete-case estimates were nearly identical. Therefore, differences among the adiposity specifications were not attributable to differential missingness or changes in sample composition.

Interpretation: Adiposity accounts for substantial overlapping variation in the HOMA-IR–hsCRP association, as indicated by attenuation of the sex-specific slopes after adiposity adjustment. However, the estimated sex difference was not eliminated by any of the three adiposity representations examined. These analyses are adjustment-sensitivity analyses and should not be interpreted as evidence of mediation, causal mechanism, or complete control of adiposity-related confounding.

Decision: Waist circumference remains the locked primary adiposity adjustment. The alternative BMI and waist-to-height specifications are retained as sensitivity analyses. The stronger interaction statistic under waist adjustment is not used as a basis for selecting waist over the alternatives; waist was already the primary specification. Stage 2A is complete and frozen.

Stage 2B — Smoking Sensitivity

Smoking status was evaluated as an expanded adjustment variable after the primary model and adiposity sensitivity analyses had been frozen. Before examining its effect on the HOMA-IR × sex interaction, smoking status was prespecified from the NHANES smoking questionnaire as:

* Never: SMQ020 = 2 (<100 cigarettes smoked in lifetime)
* Former: SMQ020 = 1 and SMQ040 = 3 (≥100 lifetime cigarettes, but not currently smoking)
* Current: SMQ020 = 1 and SMQ040 = 1 or 2 (≥100 lifetime cigarettes and currently smoking every day or some days)

Responses coded as refused, don’t know, or otherwise unclassifiable were treated as missing. No alternative smoking definitions were examined after observing the interaction results.

Smoking status was classifiable for 3,075 of the 3,079 participants in the locked primary analytic domain (99.87%). Four participants were unclassifiable (one SMQ020 = 7 and three SMQ020 = 9). To distinguish the effect of smoking adjustment from changes caused by sample composition, both the reference and smoking-adjusted models were fitted to the same smoking-complete sample of n=3,075.

Reference model:
log(hsCRP) ~ HOMA_IR * is_male + C(period) * (LBXGH + waist + AGE + Non_HDL)

Smoking-adjusted model:
log(hsCRP) ~ HOMA_IR * is_male + C(period) * (LBXGH + waist + AGE + Non_HDL) + C(smoking_status)

On the common n=3,075 sample, the reference model produced a HOMA-IR × sex interaction of F(1,40) = 9.0765, p = 0.00448. The female HOMA-IR slope was 0.06284 and the male slope was 0.00806, giving a male-minus-female slope contrast of -0.05477 (95% CI -0.09152 to -0.01803).

After adjustment for smoking status, the interaction was F(1,40) = 8.6687, p = 0.00537. The female HOMA-IR slope was 0.06321 and the male slope was 0.00985. The male-minus-female slope contrast was -0.05335 (95% CI -0.08998 to -0.01673).

Smoking adjustment therefore changed the male-minus-female contrast from -0.05477 to -0.05335, corresponding to a 2.6% reduction in its absolute magnitude. The sex-specific slopes were also essentially unchanged.

Interpretation: Adjustment for smoking status had little effect on the estimated sex difference in the HOMA-IR–hsCRP association. This analysis is interpreted as an adjustment-sensitivity analysis and not as evidence that smoking is causally unrelated to the observed association or that smoking has no biological role.

Decision: Smoking sensitivity is complete and frozen. The prespecified never/former/current definition will be retained; no alternative smoking specifications will be explored based on the observed results.

Stage 2C — Medication Use

2C. Medication use — statin sensitivity

Prespecified exposure

The principal medication sensitivity analysis evaluates current prescription
statin use. Statin use was defined from the NHANES prescription-medication
files as reported prescription medication use during the preceding 30 days
with at least one medication containing an ingredient classified in the
NHANES drug-information file as:

HMG-COA REDUCTASE INHIBITORS (STATINS)

Ingredient-level therapeutic classifications were used so that
statin-containing combination products were captured.

This definition was selected before fitting the medication-adjusted outcome
model. A broader lipid-lowering medication category was not used because it
would combine pharmacologically distinct therapies.

Medication linkage and classification QC

Prescription medication records from 2015–2016 and 2017–March 2020 were
linked to the NHANES drug-information file using RXDDRGID. Among 37,612
prescription records with a nonblank drug identifier, 37,612 (100.0%) mapped
successfully to the drug-information file. Fourteen drug identifiers
contained a statin ingredient under the prespecified therapeutic-class
definition, including both single-agent statins and statin-containing
combination products.

The locked Study #2 analytic domain was reproduced exactly (n=3,079).

Participant-level statin classification in this domain was:

* definite statin user: n=286
* definite nonuser: n=2,755
* ambiguous statin status: n=38

Ambiguous status resulted primarily from reported prescription medication
use with an unidentified medication record (e.g., 55555, 77777, or 99999);
three participants had refused/don’t-know prescription-use responses.
Participants with an identified statin were classified as statin users even
if another reported medication was unidentified.

Because an unidentified medication cannot establish absence of statin use,
participants with no detected statin but unresolved medication information
were not classified as nonusers. The 38 participants with ambiguous statin
status will therefore be excluded from the Stage 2C common-sample analysis.

This leaves n=3,041 participants with classifiable statin status.

Unweighted QC counts by sex:

* Female: 136 definite statin, 1,481 definite nonstatin, 21 ambiguous
* Male: 150 definite statin, 1,274 definite nonstatin, 17 ambiguous

These counts are descriptive QC only and are not population prevalence
estimates.

Planned Stage 2C sensitivity analysis

The medication sensitivity will compare two models fitted to the identical
common sample of n=3,041 participants.

Reference model:

log(hsCRP) ~ HOMA_IR * is_male
+ C(period) * (LBXGH + waist + AGE + Non_HDL)

Statin-adjusted model:

log(hsCRP) ~ HOMA_IR * is_male
+ C(period) * (LBXGH + waist + AGE + Non_HDL)
+ statin_use

Statin use will enter as a participant-level binary adjustment covariate.
No statin-by-sex, statin-by-HOMA-IR, or higher-order interaction will be
introduced in this sensitivity analysis.

The purpose is to determine whether adjustment for current prescription
statin use materially changes the estimated HOMA-IR × sex contrast. The
reference and statin-adjusted models will be compared on the same common
sample so that changes attributable to exclusion of participants with
ambiguous medication status can be distinguished from changes associated
with statin adjustment itself.

Status: exposure definition and QC frozen; outcome models not yet fitted.

Stage 2C results

The Stage 2C statin sensitivity analysis was conducted among the 3,041
participants with classifiable statin status. Both the reference model and
the statin-adjusted model were fitted to this identical common sample.

In the common-sample reference model, the HOMA-IR × sex interaction was:

F(1,40) = 9.261248, p = 0.004123.

The adjusted HOMA-IR slope was 0.062804 among females and 0.008032 among
males. The male-minus-female slope contrast was -0.054772
(95% CI: -0.091147 to -0.018397).

After adjustment for current prescription statin use, the HOMA-IR × sex
interaction was:

F(1,40) = 9.370413, p = 0.003929.

The adjusted HOMA-IR slope was 0.062960 among females and 0.008115 among
males. The male-minus-female slope contrast was -0.054845
(95% CI: -0.091057 to -0.018634).

Adjustment for statin use changed the male-minus-female contrast by only
-0.000074, corresponding to a 0.13% change in its absolute magnitude.
The female and male slopes were likewise essentially unchanged.

The restriction to participants with classifiable statin status itself had
negligible impact: the common-sample reference contrast (-0.054772) was
very similar to the Study #2 primary interaction-model contrast in the full
analytic sample (-0.054692).

Interpretation: adjustment for current prescription statin use did not
materially alter the estimated sex difference in the adjusted association
between HOMA-IR and hsCRP. Within the limits of this sensitivity analysis,
current statin use does not appear to account for the observed sex contrast.
This should not be interpreted as evidence that statins have no biological
effects on inflammation, lipid metabolism, or insulin resistance; rather,
accounting for measured current statin use did not materially change the
focal HOMA-IR × sex estimate in this analysis.

Status: Stage 2C complete and frozen.

2D. Renal Function

Prespecified renal specifications

Renal function was evaluated as a potential alternative explanation for the
observed sex difference in the adjusted HOMA-IR–hsCRP association.

Serum creatinine (LBXSCR) was obtained from the NHANES Standard Biochemistry
Profile. Within the locked Study #2 analytic domain (n=3,079), serum
creatinine was available for 3,075 participants (99.87%). Missingness was
minimal and balanced by sex and study period.

Renal filtration was estimated using the 2021 CKD-EPI creatinine equation.
Continuous eGFR was prespecified as the principal renal-function sensitivity.

Because the CKD-EPI equation incorporates age and sex, and sex is the focal
effect modifier in Study #2, serum creatinine was also evaluated as a
complementary directly measured renal specification. Raw serum creatinine
was strongly right-skewed in the analytic sample (skewness 15.76; kurtosis
479.71). Log transformation substantially reduced this asymmetry (skewness
0.71; kurtosis 6.09). Log serum creatinine was therefore prespecified as the
complementary renal specification before fitting any renal-adjusted outcome
model.

No participants were excluded on the basis of low eGFR or high serum
creatinine. eGFR thresholds were used for descriptive QC only.

Stage 2D results

All models were fitted to the identical renal common sample of 3,075
participants.

In the common-sample reference model, the HOMA-IR × sex interaction was:

F(1,40) = 8.913896, p = 0.004811.

The adjusted HOMA-IR slope was 0.062507 among females and 0.007442 among
males. The male-minus-female slope contrast was -0.055065
(95% CI: -0.092341 to -0.017790).

After adjustment for continuous 2021 CKD-EPI eGFR, the HOMA-IR × sex
interaction was:

F(1,40) = 9.063425, p = 0.004501.

The adjusted HOMA-IR slope was 0.062295 among females and 0.007021 among
males. The male-minus-female slope contrast was -0.055273
(95% CI: -0.092380 to -0.018167).

eGFR adjustment changed the contrast by -0.000208, corresponding to only
a 0.38% change in its absolute magnitude.

In the complementary log-serum-creatinine model, the HOMA-IR × sex
interaction was:

F(1,40) = 9.049677, p = 0.004529.

The adjusted HOMA-IR slope was 0.062329 among females and 0.007107 among
males. The male-minus-female slope contrast was -0.055222
(95% CI: -0.092323 to -0.018122).

Log-creatinine adjustment changed the contrast by -0.000157, corresponding
to only a 0.29% change in its absolute magnitude.

The restriction to participants with available renal data also had little
impact: the common-sample reference contrast (-0.055065) was very similar
to the Study #2 primary interaction-model contrast in the full analytic
sample (-0.054692).

Interpretation: adjustment for renal function did not materially alter the
estimated sex difference in the adjusted association between HOMA-IR and
hsCRP. Results were nearly identical when renal status was represented by
continuous eGFR or by log-transformed serum creatinine. Within the limits
of these sensitivity analyses, measured renal function does not appear to
account for the observed sex contrast. This should not be interpreted as
evidence that renal physiology has no biological relationship with
inflammation or insulin resistance.

Status: Stage 2D complete and frozen.

Stage 2E — Liver / Metabolic Markers

2E. Hepatic/Metabolic Biochemistry

Prespecified specification and QC

Hepatic/metabolic biochemistry was evaluated as a potential alternative
explanation for the observed sex difference in the adjusted HOMA-IR–hsCRP
association.

ALT (LBXSATSI), AST (LBXSASSI), and GGT (LBXSGTSI) were examined during
predictor-only QC before fitting any Stage 2E outcome model. Coverage within
the locked Study #2 analytic domain was 99.84% for ALT, 99.64% for AST,
and 99.84% for GGT.

ALT was selected prospectively as the Stage 2E hepatic/metabolic sensitivity.
Raw ALT was strongly right-skewed in the analytic sample (skewness 14.90;
kurtosis 464.73). Natural-log transformation substantially reduced this
asymmetry (skewness 0.70; kurtosis 1.38). Log-transformed ALT was therefore
frozen as the Stage 2E specification before fitting the hsCRP outcome model.

AST and GGT were examined during predictor QC but were not added as
additional outcome-model sensitivities. This decision was made before
observing the Stage 2E outcome results in order to retain a targeted
sensitivity analysis rather than proliferating related hepatic biomarkers.

No participants were excluded on the basis of high or low ALT values.

Stage 2E results

Both models were fitted to the identical ALT common sample of 3,074
participants.

In the common-sample reference model, the HOMA-IR × sex interaction was:

F(1,40) = 8.910135, p = 0.004819.

The adjusted HOMA-IR slope was 0.062520 among females and 0.007466 among
males. The male-minus-female slope contrast was -0.055054
(95% CI: -0.092331 to -0.017778).

After adjustment for log-transformed ALT, the HOMA-IR × sex interaction was:

F(1,40) = 8.809834, p = 0.005041.

The adjusted HOMA-IR slope was 0.063278 among females and 0.008413 among
males. The male-minus-female slope contrast was -0.054866
(95% CI: -0.092225 to -0.017506).

Log-ALT adjustment changed the male-minus-female contrast by +0.000189,
corresponding to a 0.34% reduction in its absolute magnitude.

The restriction to participants with available ALT also had negligible
impact: the common-sample reference contrast (-0.055054) was very similar
to the Study #2 primary interaction-model contrast in the full analytic
sample (-0.054692).

Interpretation: adjustment for log-transformed ALT did not materially alter
the estimated sex difference in the adjusted association between HOMA-IR
and hsCRP. Within the limits of this sensitivity analysis, variation in ALT
does not appear to account for the observed sex contrast. This should not
be interpreted as evidence that hepatic physiology or metabolic liver
processes have no biological relationship with insulin resistance or
inflammation.

Status: Stage 2E complete and frozen.

### 2F. Current Pregnancy

#### Rationale and prespecified sensitivity definition

Because pregnancy produces substantial temporary changes in insulin sensitivity,
endocrine physiology, body composition, lipid metabolism, and inflammatory
physiology, currently pregnant participants represent a physiologically
distinct state that could potentially influence the estimated sex difference
in the HOMA-IR–hsCRP association.

Current pregnancy status was identified using the NHANES examination pregnancy
status variable RIDEXPRG from the demographic files. Pregnancy status was
evaluated only as an exclusion sensitivity and was not added to the primary
model as an adjustment covariate.

Before fitting any pregnancy-sensitivity outcome model, pregnancy status was
examined within the locked Study #2 analytic domain (n=3,079).

Among 1,638 female participants:

- 857 were classified as not pregnant;
- 48 were classified as currently pregnant;
- 15 had pregnancy status that could not be determined; and
- 718 had no reported RIDEXPRG value.

All 718 females without a reported RIDEXPRG value were age 45 years or older.
Among the 920 females ages 20–44, pregnancy status was therefore classified
for every participant: 857 not pregnant, 48 pregnant, and 15 cannot determine.

The principal Stage 2F sensitivity was prespecified as exclusion of participants
with confirmed current pregnancy (RIDEXPRG=1), while retaining participants
whose pregnancy status could not be determined. This produces an expected
analysis sample of n=3,031.

A complementary conservative sensitivity was prespecified excluding both
confirmed pregnancy (RIDEXPRG=1) and cannot-determine pregnancy status
(RIDEXPRG=3), producing an expected analysis sample of n=3,016.

These definitions were frozen before examining the HOMA-IR × sex interaction
under either pregnancy exclusion.

No exclusion based on prior pregnancy, parity, history of gestational diabetes,
or postpartum status was included in this stage. Such variables represent
distinct reproductive-history questions and do not provide longitudinal
before-versus-after pregnancy comparisons in the cross-sectional NHANES design.

#### Stage 2F results

The principal pregnancy sensitivity excluded the 48 participants with confirmed
current pregnancy (RIDEXPRG=1), leaving n=3,031 participants.

After exclusion of confirmed current pregnancy, the HOMA-IR × sex interaction
was:

F(1,40) = 9.044473, p = 0.004539.

The adjusted HOMA-IR slope was 0.066633 among females and 0.009884 among males,
giving a male-minus-female slope contrast of -0.056748.

Relative to the full Study #2 primary interaction-model contrast of -0.054692,
the absolute magnitude of the sex contrast increased modestly (approximately
3.8%). Exclusion of confirmed pregnancy therefore did not materially attenuate
the estimated sex difference.

In the complementary conservative sensitivity excluding both confirmed
pregnancy (RIDEXPRG=1) and cannot-determine pregnancy status (RIDEXPRG=3),
n=3,016 participants remained. The HOMA-IR × sex interaction was:

F(1,40) = 8.943839, p = 0.004747.

The adjusted HOMA-IR slope was 0.066827 among females and 0.010080 among males,
giving a male-minus-female slope contrast of -0.056747.

The nearly identical contrasts under the two pregnancy exclusions indicate that
treatment of the 15 participants with cannot-determine pregnancy status had
negligible influence on the result.

#### Interpretation

Exclusion of participants with confirmed current pregnancy did not materially
alter the estimated sex difference in the adjusted HOMA-IR–hsCRP association.
The same conclusion was obtained under the more conservative exclusion of both
confirmed pregnancy and cannot-determine pregnancy status.

These analyses address current pregnancy only. They should not be interpreted
as accounting for menopausal status, prior pregnancy, parity, history of
gestational diabetes, exogenous hormone use, or other reproductive or endocrine
states.

#### Decision

Stage 2F is complete and frozen. Confirmed current pregnancy will not be
excluded from the locked Study #2 primary analytic population. The pregnancy
exclusions are retained as sensitivity analyses, and no additional
pregnancy-related exclusions will be pursued for Study #2.

### Stage 2G — Alcohol Adjustment Sensitivity

#### Rationale

Alcohol use was evaluated as an additional sensitivity covariate because
alcohol exposure may be associated with systemic inflammatory markers and
may differ by sex. Alcohol was not added to the Study #2 primary interaction
model. The purpose of this analysis was specifically to determine whether
additional adjustment for broad alcohol-use status materially changed the
estimated HOMA-IR × sex interaction.

The Study #2 primary interaction model remained:

    log_hs_CRP ~ HOMA_IR * is_male
                 + C(period) * (LBXGH + waist + AGE + Non_HDL)

The alcohol sensitivity model added:

    + C(alcohol_status)

Both the reference and alcohol-adjusted models were fit on the identical
alcohol-classifiable sample so that changes in the focal interaction estimate
could be attributed to alcohol adjustment rather than to differences in
sample composition.

#### Exposure reconstruction and pre-outcome specification

Alcohol questionnaire structure differed between the 2015–2016 and
2017–March 2020 cycles. Alcohol exposure was therefore reconstructed from
cycle-specific questionnaire variables before examining HOMA-IR/hsCRP
results.

SAS/XPORT pseudo-zero values (~5.397605e-79) were normalized to numeric zero
before classification.

The categorical alcohol specification was frozen before fitting any
HOMA-IR/hsCRP alcohol-sensitivity model.

Final categories were:

- Zero: documented zero past-year alcohol exposure.
- Very-low lifetime: 2015–2016 participants reporting fewer than 12 drinks
  in any year and fewer than 12 drinks in their lifetime, but for whom exact
  past-year volume could not be reconstructed.
- Positive: documented positive past-year drinking frequency.
- Unknown/unclassifiable: insufficient information to assign one of the
  above categories.

Locked-domain counts were:

- Zero: 555
- Very-low lifetime: 163
- Positive: 2,201
- Unknown/unclassifiable: 160
- Total: 3,079

The alcohol-classifiable sensitivity sample was therefore n=2,919.

During pre-outcome QC, two locked-cohort participants in the later survey
period were found to have valid positive past-year drinking-frequency
responses (ALQ121 > 0) but ALQ130=999, indicating unavailable/invalid typical
drinks per drinking day. Earlier volume-oriented QC had classified these two
participants as unknown because exact drinks/week could not be estimated.
For the subsequently specified categorical alcohol-status sensitivity, they
were classified as Positive because past-year drinking status was directly
observed even though drinking volume was not. This refinement was made before
fitting any HOMA-IR/hsCRP alcohol-sensitivity model.

No continuous drinks/week term, dose categories, polynomial terms, splines,
or alternative alcohol cutpoints were selected after examining the outcome.

#### Results

Both models used the same n=2,919 participants and retained design
denominator df=40.

Same-sample reference model:

- HOMA-IR × sex: F(1,40)=12.381615, p=0.001097
- Female HOMA-IR slope: 0.075707
- Male HOMA-IR slope: 0.008394
- Male-minus-female slope difference: -0.067313
- 95% CI for slope difference: [-0.105976, -0.028650]

Alcohol-adjusted model:

- HOMA-IR × sex: F(1,40)=12.402796, p=0.001088
- Female HOMA-IR slope: 0.075024
- Male HOMA-IR slope: 0.008492
- Male-minus-female slope difference: -0.066532
- 95% CI for slope difference: [-0.104714, -0.028350]

Change after alcohol adjustment:

- Change in male-minus-female slope contrast:
  +0.000781
- Change in absolute contrast magnitude:
  -1.16%

#### Interpretation

Additional adjustment for broad categorical alcohol-use status had little
effect on the estimated sex difference in the HOMA-IR–hsCRP association.

The male-minus-female slope contrast changed from -0.067313 in the
same-sample reference model to -0.066532 after alcohol adjustment, a 1.16%
reduction in absolute magnitude. Female and male HOMA-IR slopes were also
nearly unchanged.

The difference between the full-cohort primary interaction estimate and the
same-sample reference estimate in this sensitivity should not be attributed
to alcohol adjustment. The same-sample reference model already reflects
restriction to participants with classifiable alcohol status. The relevant
assessment of alcohol adjustment is therefore the within-sample comparison
between the reference and alcohol-adjusted models.

These results do not establish that alcohol has no biological relationship
with insulin resistance or systemic inflammation. Rather, under the
prespecified categorical sensitivity specification, additional adjustment
for alcohol-use status did not materially account for the observed
HOMA-IR × sex interaction.

#### Decision

Stage 2G is complete and frozen.

Alcohol will not be added to the Study #2 primary interaction model.

No additional alcohol specifications will be pursued for Study #2 unless a
specific methodological concern subsequently requires them.

Script:

    scripts/sex_alcohol_sensitivity.py

Machine-readable output:

    output/sex_alcohol_sensitivity/alcohol_sensitivity.csv

Console output:

    output/sex_alcohol_sensitivity.txt

### 2H hsCRP ≤10 mg/L

**Purpose.** Evaluate whether the Study #2 primary sex interaction is sensitive
to the upper range of the hsCRP distribution. The 10 mg/L threshold was not
selected from the Study #2 results; it had already been used as a sensitivity
restriction in Study #1.

**Pre-outcome QC.**
- Locked Study #2 primary cohort: n = 3,079.
- hsCRP was available and positive for all 3,079 participants.
- hsCRP >10 mg/L: n = 184.
- hsCRP <=10 mg/L: n = 2,895.
- The restriction was frozen at hsCRP <=10 mg/L before fitting the sensitivity
  model.
- No alternate hsCRP thresholds were explored.

**Primary-model reference (full cohort, n = 3,079; design df = 40):**
- HOMA-IR × sex: F(1,40) = 9.053325, p = 0.004522.
- Female HOMA-IR slope = 0.062864.
- Male HOMA-IR slope = 0.008172.
- Male-minus-female slope difference = -0.054692.
- 95% CI for male-minus-female difference: [-0.091428, -0.017955].

**hsCRP <=10 mg/L sensitivity (n = 2,895; design df = 40):**
- HOMA-IR × sex: F(1,40) = 3.236757, p = 0.079551.
- Female HOMA-IR slope = 0.043762.
- Male HOMA-IR slope = 0.016839.
- Male-minus-female slope difference = -0.026922.
- 95% CI for male-minus-female difference: [-0.057167, 0.003322].

**Change from the full primary model.**
- Male-minus-female interaction contrast changed from -0.054692 to -0.026922.
- Absolute interaction magnitude was attenuated by approximately 50.8%.
- The female slope decreased from 0.062864 to 0.043762.
- The male slope increased from 0.008172 to 0.016839.
- Thus, the attenuation reflected convergence of the two sex-specific slopes,
  rather than a change confined to only one sex.

**Interpretation.**
The estimated sex difference was materially attenuated after restriction to
hsCRP <=10 mg/L. This is not merely a change in whether the interaction crossed
a conventional significance threshold: the estimated interaction contrast
itself decreased by approximately one-half. The Study #2 primary interaction is
therefore partly sensitive to the upper range of the hsCRP distribution.

The hsCRP >10 mg/L group should not be interpreted automatically as representing
acute infection. hsCRP is a nonspecific marker of systemic inflammation, and
values >10 mg/L may occur during acute inflammatory states but may also reflect
persistent or chronic inflammatory processes. The sensitivity restriction is
therefore treated as a robustness analysis rather than as a validated exclusion
of participants with acute infection.

#### Descriptive diagnostic of participants with hsCRP >10 mg/L

Because the hsCRP <=10 mg/L restriction materially attenuated the interaction,
a descriptive diagnostic was performed to characterize the 184 participants
with hsCRP >10 mg/L. This diagnostic was not used to define new exclusions,
select another hsCRP threshold, or construct an alternative primary model.

CBC data were available for all 184 participants from the corresponding NHANES
cycles. A small CBC panel was frozen before examining its values:
- total white blood cell count (WBC),
- absolute neutrophil count,
- neutrophil percentage,
- absolute lymphocyte count,
- absolute monocyte count,
- derived neutrophil-to-lymphocyte ratio (NLR).

Among the 184 participants with hsCRP >10 mg/L:
- 128 were female and 56 were male.
- The unweighted prevalence of hsCRP >10 mg/L was 7.81% among females and
  3.89% among males.
- Survey-weighted prevalence was 5.40% among females and 4.24% among males.

CBC measures showed a coherent descriptive shift toward greater immune
activation in the hsCRP >10 mg/L group in both sexes.

Female medians, hsCRP <=10 versus >10 mg/L:
- WBC: 6.2 vs 7.5.
- Absolute neutrophils: 3.5 vs 4.45.
- Neutrophil percentage: 57.4% vs 60.0%.
- NLR: 1.82 vs 2.13.

Male medians, hsCRP <=10 versus >10 mg/L:
- WBC: 6.3 vs 7.2.
- Absolute neutrophils: 3.5 vs 4.5.
- Neutrophil percentage: 56.6% vs 61.4%.
- NLR: 1.84 vs 2.39.

These CBC patterns provide corroborating evidence that participants with
hsCRP >10 mg/L had greater contemporaneous immune/inflammatory activation on
average. They do not establish infection, and the typical CBC values did not
justify classifying the entire hsCRP >10 mg/L group as having an acute
infectious process.

#### HOMA-IR distribution and interaction geometry

The HOMA-IR distribution differed notably by hsCRP group among females but not
in the same way among males.

Female HOMA-IR:
- hsCRP <=10 mg/L: median 1.905, Q75 2.871, P90 4.510.
- hsCRP >10 mg/L: median 3.005, Q75 4.666, P90 7.403.

Male HOMA-IR:
- hsCRP <=10 mg/L: median 2.051, Q75 3.354, P90 5.128.
- hsCRP >10 mg/L: median 2.044, Q75 4.125, P90 6.767.

Scatterplots of raw HOMA-IR against log(hsCRP), displayed separately by sex,
showed that female observations with hsCRP >10 mg/L extended farther into the
higher-HOMA portion of the joint distribution. In contrast, high-hsCRP male
observations were more concentrated at lower HOMA-IR values, while numerous
high-HOMA male observations remained below the hsCRP = 10 mg/L threshold.

This geometry provides a descriptive statistical explanation for the Stage 2H
attenuation. Restricting hsCRP to <=10 mg/L preferentially removes observations
from the upper-right portion of the female HOMA-IR–hsCRP distribution, which
would be expected to flatten the female slope. At the same time, removal of
high-hsCRP male observations without a comparable rightward HOMA-IR shift can
allow the male slope to increase. This is consistent with the observed
convergence of the sex-specific slopes.

Two diagnostic scatterplots were retained:
1. a full-range plot displaying the complete HOMA-IR distribution; and
2. a zoomed HOMA-IR 0–12 view for visualization of approximately the central
   99% of the HOMA-IR distribution.

The zoom was a display restriction only; no observations were removed from any
analysis.

**Decision.**
Stage 2H is frozen. The hsCRP <=10 mg/L result will be reported transparently
as an important sensitivity of the primary interaction. No additional hsCRP
cutoffs, CBC-based exclusions, infection classifications, or alternative
outcome models will be pursued to recover or strengthen the interaction.
The diagnostic supports describing the interaction as partly dependent on the
upper hsCRP range, while leaving unresolved whether those high-inflammatory
states reflect acute inflammation, persistent inflammation, metabolic
inflammation, or a mixture of these processes.

## Deferred Reproductive-State / Menopause Exploration

### Rationale

Because the Study #2 primary analysis identified a sex difference in the
adjusted HOMA-IR–hsCRP association, reproductive and menopausal physiology was
considered as a potentially important biological direction for further
investigation.

This exploration was undertaken to determine whether the available NHANES
reproductive-health variables could support a sufficiently defensible
reproductive-state classification before examining HOMA-IR–hsCRP outcome
associations by reproductive state.

### Predictor-only QC

Reproductive-health questionnaire variables were examined among the 1,638
female participants in the locked Study #2 analytic cohort.

Variables considered included:

- recent menstrual-period status (RHQ031);
- hysterectomy history (RHD280); and
- bilateral oophorectomy status (RHQ305).

Hysterectomy was not treated as equivalent to bilateral oophorectomy because
the two procedures have different implications for ovarian function and
endocrine physiology.

Among female participants:

- RHQ031 indicated recent menses in 900 participants and no recent menses in
  640, with 98 missing or otherwise unavailable;
- hysterectomy was reported by 241 participants and not reported by 1,288,
  with remaining responses missing or otherwise unresolved; and
- bilateral oophorectomy was reported by 120 participants and not reported by
  1,399, with remaining responses missing or otherwise unresolved.

Among women reporting no recent menses, no hysterectomy, and no bilateral
oophorectomy, substantial age heterogeneity remained. This group included
women younger than 45 years, including women younger than 40 years. Therefore,
absence of recent menses alone could not be assumed to represent natural
menopause.

A possible high-specificity reproductive-state framework was considered:

- Cycling: recent menses reported, no hysterectomy, and no bilateral
  oophorectomy.
- Natural postmenopausal phenotype: no recent menses, no hysterectomy, no
  bilateral oophorectomy, and age >=45 years.
- Surgical menopause phenotype: bilateral oophorectomy reported.
- Other/indeterminate: all remaining combinations.

Current pregnancy would require separate handling in any future
reproductive-state analysis.

### Decision

No HOMA-IR–hsCRP outcome model was fitted using these reproductive-state
categories for Study #2.

Menopausal/reproductive state was not part of the locked Study #2 primary
question, and the available cross-sectional questionnaire information does not
provide a uniquely determined menopausal classification for all women.
Constructing and comparing multiple reproductive-state definitions after
observing the Study #2 sex interaction would also risk turning a biologically
motivated follow-up question into a post-hoc model search.

Reproductive and menopausal physiology is therefore deferred as a distinct
future research question rather than added as another Study #2 sensitivity.
The Study #2 Discussion may identify reproductive/endocrine physiology as a
plausible biological direction that was not resolved by the present analysis.

No reproductive-state outcome results exist for Study #2, and no claim should
be made that menopause explains the observed sex difference.

Status: exploratory exposure/QC work documented; reproductive-state outcome
analysis deliberately deferred.
