#!/usr/bin/env Rscript

# Cross-software validation of Study #2 primary NHANES sex-interaction analysis
# -----------------------------------------------------------------------------
# Purpose:
#   Independently reproduce the locked Study #2 primary Python result using
#   R's established 'survey' package.
#
#   This is a validation implementation, not a translation of the Python
#   survey-estimation code. It starts from the raw NHANES XPT files and uses
#   survey::svydesign(), subset(), svyglm(), and design-based contrasts.
#
# Study #2 primary model:
#
#   log_hs_CRP ~ HOMA_IR * is_male
#                + period * (LBXGH + waist + AGE + Non_HDL)
#
# Run:
#   Rscript validation/validate_primary_sex_R.R \
#     --data-dir /path/to/data_raw
#
# Required packages:
#   haven, survey
#
# Install once if needed:
#   install.packages(c("haven", "survey"))

options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

required <- c("haven", "survey")
missing_pkgs <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall once with:\n  install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(haven)
  library(survey)
})

# ------------------------------------------------------------
# Command-line argument parsing
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2 || args[1] != "--data-dir") {
  stop(
    "Usage:\n",
    "  Rscript validation/validate_primary_sex_R.R ",
    "--data-dir /path/to/data_raw",
    call. = FALSE
  )
}

data_dir <- normalizePath(args[2], mustWork = TRUE)

# ------------------------------------------------------------
# Locked study specification
# ------------------------------------------------------------

MIN_AGE <- 20
HBA1C_MAX <- 5.7
WEIGHT_TOL <- 1e-12
TOTAL_YEARS <- 5.2

# Frozen Python targets from Study #2.
#
# These values are used ONLY after R independently constructs the cohort,
# survey design, and fitted model.

expected <- list(
  n = 3079L,
  design_df = 40,

  interaction_F = 9.053325,
  interaction_p = 0.004521562,

  female_beta = 0.062864,
  female_se = 0.018001,
  female_ci_low = 0.026483,
  female_ci_high = 0.099245,
  female_p = 0.001184,

  male_beta = 0.008172,
  male_se = 0.011743,
  male_ci_low = -0.015561,
  male_ci_high = 0.031905,
  male_p = 0.490492,

  male_minus_female = -0.054692,
  difference_ci_low = -0.091428,
  difference_ci_high = -0.017955
)

# Recorded Python values are rounded to approximately six decimal places.
# This tolerance permits only harmless differences at that recorded precision.
TOL <- 1e-6

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

read_xpt_checked <- function(path) {

  if (!file.exists(path)) {
    stop(
      "Required NHANES file not found: ",
      path,
      call. = FALSE
    )
  }

  d <- haven::read_xpt(path)

  if (!("SEQN" %in% names(d))) {
    stop(
      basename(path),
      " is missing SEQN",
      call. = FALSE
    )
  }

  if (anyDuplicated(d$SEQN)) {
    stop(
      basename(path),
      " contains duplicate SEQN values",
      call. = FALSE
    )
  }

  d
}


nhanes_filename <- function(key, prefix, suffix) {

  stems <- c(
    demo  = "DEMO",
    bmx   = "BMX",
    diq   = "DIQ",
    fast  = "FASTQX",
    ghb   = "GHB",
    glu   = "GLU",
    ins   = "INS",
    hscrp = "HSCRP",
    tchol = "TCHOL",
    hdl   = "HDL"
  )

  paste0(prefix, stems[[key]], suffix)
}


left_merge_one_to_one <- function(x, y, cols) {

  keep <- intersect(cols, names(y))
  y2 <- y[, keep, drop = FALSE]

  if (anyDuplicated(y2$SEQN)) {
    stop(
      "Duplicate SEQN in merge source",
      call. = FALSE
    )
  }

  # Preserve DEMO row order exactly.
  idx <- match(x$SEQN, y2$SEQN)

  for (nm in setdiff(names(y2), "SEQN")) {
    x[[nm]] <- y2[[nm]][idx]
  }

  x
}


build_period <- function(
  label,
  prefix,
  suffix,
  fasting_weight,
  duration_years
) {

  keys <- c(
    "demo",
    "bmx",
    "diq",
    "fast",
    "ghb",
    "glu",
    "ins",
    "hscrp",
    "tchol",
    "hdl"
  )

  files <- setNames(
    lapply(
      keys,
      function(k) {
        read_xpt_checked(
          file.path(
            data_dir,
            nhanes_filename(k, prefix, suffix)
          )
        )
      }
    ),
    keys
  )

  demo_cols <- c(
    "SEQN",
    "RIDAGEYR",
    "RIAGENDR",
    "SDMVSTRA",
    "SDMVPSU"
  )

  d <- as.data.frame(
    files$demo[, demo_cols]
  )

  selections <- list(

    bmx = c(
      "SEQN",
      "BMXWAIST",
      "BMXHT",
      "BMXBMI",
      "BMXHIP"
    ),

    diq = c(
      "SEQN",
      "DIQ010"
    ),

    fast = c(
      "SEQN",
      "PHAFSTHR",
      "PHAFSTMN"
    ),

    ghb = c(
      "SEQN",
      "LBXGH"
    ),

    glu = c(
      "SEQN",
      fasting_weight,
      "LBXGLU"
    ),

    ins = c(
      "SEQN",
      "LBXIN"
    ),

    hscrp = c(
      "SEQN",
      "LBXHSCRP"
    ),

    tchol = c(
      "SEQN",
      "LBXTC"
    ),

    hdl = c(
      "SEQN",
      "LBDHDD"
    )
  )

  for (k in names(selections)) {

    d <- left_merge_one_to_one(
      d,
      as.data.frame(files[[k]]),
      selections[[k]]
    )
  }

  # ----------------------------------------------------------
  # Harmonized variables
  # ----------------------------------------------------------

  d$period <- label

  d$raw_fasting_weight <- d[[fasting_weight]]

  d$analysis_weight <-
    d$raw_fasting_weight *
    (duration_years / TOTAL_YEARS)

  d$AGE <- d$RIDAGEYR

  # NHANES RIAGENDR:
  #   1 = Male
  #   2 = Female
  #
  # Study #2 coding:
  #   female = 0
  #   male   = 1

  d$is_male <- ifelse(
    d$RIAGENDR == 1,
    1,
    ifelse(
      d$RIAGENDR == 2,
      0,
      NA_real_
    )
  )

  d$waist <- d$BMXWAIST
  d$height <- d$BMXHT
  d$BMI <- d$BMXBMI

  if ("BMXHIP" %in% names(d)) {
    d$WHR <- d$BMXWAIST / d$BMXHIP
  } else {
    d$WHR <- NA_real_
  }

  d$waist_to_height <-
    d$waist / d$height

  d$hs_CRP <- d$LBXHSCRP

  d$log_hs_CRP <- ifelse(
    d$hs_CRP > 0,
    log(d$hs_CRP),
    NA_real_
  )

  d$HOMA_IR <-
    (d$LBXIN * d$LBXGLU) / 405

  d$Non_HDL <-
    d$LBXTC - d$LBDHDD

  d$fasting_hours <-
    d$PHAFSTHR +
    d$PHAFSTMN / 60

  # Unique masked design IDs across periods.
  d$stratum_u <-
    paste(
      d$period,
      d$SDMVSTRA,
      sep = ":"
    )

  d$psu_u <-
    paste(
      d$period,
      d$SDMVPSU,
      sep = ":"
    )

  d
}


# ------------------------------------------------------------
# Header
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("R CROSS-SOFTWARE VALIDATION\n")
cat("STUDY #2 PRIMARY HOMA-IR x SEX INTERACTION\n")
cat("============================================================\n\n")

cat(
  "R version:      ",
  R.version.string,
  "\n",
  sep = ""
)

cat(
  "survey version: ",
  as.character(packageVersion("survey")),
  "\n",
  sep = ""
)

cat(
  "haven version:  ",
  as.character(packageVersion("haven")),
  "\n\n",
  sep = ""
)

# ------------------------------------------------------------
# Build data independently from raw NHANES XPT files
# ------------------------------------------------------------

d15 <- build_period(
  label = "2015-16",
  prefix = "",
  suffix = "_I.xpt",
  fasting_weight = "WTSAF2YR",
  duration_years = 2.0
)

d20 <- build_period(
  label = "2017-Mar2020",
  prefix = "P_",
  suffix = ".xpt",
  fasting_weight = "WTSAFPRP",
  duration_years = 3.2
)

# Harmonize period-specific columns before stacking.
# Some NHANES variables are available only in one period.
# They are not required by the primary Study #2 model, but
# base R rbind() requires identical columns.

all_cols <- union(
  names(d15),
  names(d20)
)


add_missing_columns <- function(x, all_cols) {

  missing_cols <-
    setdiff(
      all_cols,
      names(x)
    )

  for (nm in missing_cols) {
    x[[nm]] <- NA
  }

  x[, all_cols, drop = FALSE]
}


d15 <- add_missing_columns(
  d15,
  all_cols
)

d20 <- add_missing_columns(
  d20,
  all_cols
)

d <- rbind(
  d15,
  d20
)

# Explicit period reference level.
d$period <- factor(
  d$period,
  levels = c(
    "2015-16",
    "2017-Mar2020"
  )
)

# ------------------------------------------------------------
# Eligibility / participant flow
# ------------------------------------------------------------

age_ok <-
  !is.na(d$AGE) &
  d$AGE >= MIN_AGE

diab_ok <-
  !is.na(d$DIQ010) &
  d$DIQ010 == 2

hba1c_ok <-
  !is.na(d$LBXGH) &
  d$LBXGH < HBA1C_MAX

fasting_w_ok <-
  !is.na(d$raw_fasting_weight) &
  d$raw_fasting_weight > WEIGHT_TOL

labs_ok <-
  !is.na(d$LBXGLU) &
  !is.na(d$LBXIN)

positive_ok <-
  !is.na(d$HOMA_IR) &
  !is.na(d$hs_CRP) &
  d$HOMA_IR > 0 &
  d$hs_CRP > 0

covariate_complete <-
  complete.cases(
    d[, c(
      "LBXGH",
      "waist",
      "AGE",
      "is_male",
      "Non_HDL"
    )]
  )

m0 <- rep(
  TRUE,
  nrow(d)
)

m1 <- m0 & age_ok
m2 <- m1 & diab_ok
m3 <- m2 & hba1c_ok
m4 <- m3 & fasting_w_ok
m5 <- m4 & labs_ok
m6 <- m5 & positive_ok
m7 <- m6 & covariate_complete

flow <- data.frame(

  stage = c(
    "NHANES 2015-March 2020 participants",
    "Age >= 20",
    "DIQ010 == 2",
    "HbA1c < 5.7%",
    "Qualifying fasting-subsample weight",
    "Valid fasting glucose and insulin",
    "Positive HOMA-IR and hsCRP",
    "Complete primary covariates"
  ),

  remaining_n = c(
    sum(m0),
    sum(m1),
    sum(m2),
    sum(m3),
    sum(m4),
    sum(m5),
    sum(m6),
    sum(m7)
  )
)

flow$excluded_at_stage <-
  c(
    NA,
    -diff(flow$remaining_n)
  )

cat("PARTICIPANT FLOW\n")
print(
  flow,
  row.names = FALSE
)
cat("\n")

d$primary_domain <- m7

# ------------------------------------------------------------
# Survey design
# ------------------------------------------------------------

# Construct the design on the qualifying fasting-subsample
# frame, then estimate the final analytic population as a
# domain within that design.
#
# This preserves fasting-subsample participants outside the
# final analytic domain for domain variance estimation.

fasting_design_frame <- d[
  !is.na(d$raw_fasting_weight) &
    d$raw_fasting_weight > WEIGHT_TOL &
    !is.na(d$analysis_weight),
  ,
  drop = FALSE
]

fasting_design <- survey::svydesign(
  ids = ~psu_u,
  strata = ~stratum_u,
  weights = ~analysis_weight,
  data = fasting_design_frame,
  nest = TRUE
)

primary_design <-
  subset(
    fasting_design,
    primary_domain
  )

analytic_n <-
  sum(d$primary_domain)

design_df <-
  survey::degf(primary_design)

cat("DESIGN CHECKS\n")

cat(
  "Analytic n: ",
  analytic_n,
  "\n",
  sep = ""
)

cat(
  "Survey df:  ",
  design_df,
  "\n\n",
  sep = ""
)

# ------------------------------------------------------------
# Locked Study #2 primary interaction model
# ------------------------------------------------------------

# Python/Patsy specification:
#
# log_hs_CRP ~
#   HOMA_IR * is_male +
#   C(period) * (LBXGH + waist + AGE + Non_HDL)

primary_formula <- log_hs_CRP ~
  HOMA_IR * is_male +
  period * (
    LBXGH +
    waist +
    AGE +
    Non_HDL
  )

fit <- survey::svyglm(
  primary_formula,
  design = primary_design,
  family = gaussian()
)

cat("PRIMARY MODEL COEFFICIENTS\n")
print(
  coef(summary(fit))
)
cat("\n")

# ------------------------------------------------------------
# Locate focal coefficient names defensively
# ------------------------------------------------------------

b <- coef(fit)
V <- vcov(fit)

needed_names <- c(
  "HOMA_IR",
  "HOMA_IR:is_male"
)

missing_names <-
  setdiff(
    needed_names,
    names(b)
  )

if (length(missing_names) > 0) {

  stop(
    "Expected coefficient name(s) not found: ",
    paste(
      missing_names,
      collapse = ", "
    ),
    "\nAvailable names:\n",
    paste(
      names(b),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Formal HOMA-IR x sex interaction
# ------------------------------------------------------------

# Because the interaction contains one parameter, the Wald
# statistic is t^2 and is evaluated as F(1, design_df).

interaction_beta <-
  unname(
    b["HOMA_IR:is_male"]
  )

interaction_se <-
  sqrt(
    V[
      "HOMA_IR:is_male",
      "HOMA_IR:is_male"
    ]
  )

interaction_F <-
  (interaction_beta /
     interaction_se)^2

interaction_p <-
  pf(
    interaction_F,
    df1 = 1,
    df2 = design_df,
    lower.tail = FALSE
  )

cat("FORMAL HOMA-IR x SEX INTERACTION\n")

cat(
  sprintf(
    "F(1,%d) = %.12f, p = %.12f\n\n",
    design_df,
    interaction_F,
    interaction_p
  )
)

# ------------------------------------------------------------
# Linear-contrast helper
# ------------------------------------------------------------

contrast_stats <- function(
  L,
  label
) {

  est <-
    sum(
      L * b
    )

  se <-
    sqrt(
      as.numeric(
        t(L) %*%
          V %*%
          L
      )
    )

  tval <-
    est / se

  pval <-
    2 *
    pt(
      abs(tval),
      df = design_df,
      lower.tail = FALSE
    )

  crit <-
    qt(
      0.975,
      df = design_df
    )

  data.frame(
    quantity = label,
    beta = est,
    SE = se,
    CI_low = est - crit * se,
    CI_high = est + crit * se,
    t = tval,
    df = design_df,
    p = pval,
    row.names = NULL
  )
}

# ------------------------------------------------------------
# Female slope
# ------------------------------------------------------------

# Female is the reference category:
#
# beta_female = beta_HOMA

L_female <-
  setNames(
    rep(
      0,
      length(b)
    ),
    names(b)
  )

L_female["HOMA_IR"] <- 1

# ------------------------------------------------------------
# Male slope
# ------------------------------------------------------------

# beta_male =
#   beta_HOMA + beta_HOMA:male

L_male <- L_female

L_male[
  "HOMA_IR:is_male"
] <- 1

# ------------------------------------------------------------
# Male-minus-female slope difference
# ------------------------------------------------------------

# beta_male - beta_female =
#   beta_HOMA:male

L_difference <-
  setNames(
    rep(
      0,
      length(b)
    ),
    names(b)
  )

L_difference[
  "HOMA_IR:is_male"
] <- 1

female <-
  contrast_stats(
    L_female,
    "Female HOMA-IR slope"
  )

male <-
  contrast_stats(
    L_male,
    "Male HOMA-IR slope"
  )

difference <-
  contrast_stats(
    L_difference,
    "Male-minus-female slope difference"
  )

contrasts <-
  rbind(
    female,
    male,
    difference
  )

cat("SEX-SPECIFIC HOMA-IR SLOPES AND CONTRAST\n")

print(
  contrasts,
  digits = 12,
  row.names = FALSE
)

cat("\n")

# ------------------------------------------------------------
# Validation against frozen Python results
# ------------------------------------------------------------

check_num <- function(
  label,
  expected_value,
  obtained_value,
  tol = TOL
) {

  diff <-
    abs(
      expected_value -
      obtained_value
    )

  pass <-
    is.finite(diff) &&
    diff <= tol

  data.frame(
    quantity = label,
    expected = expected_value,
    R_obtained = obtained_value,
    abs_difference = diff,
    tolerance = tol,
    status = ifelse(
      pass,
      "PASS",
      "FAIL"
    ),
    row.names = NULL
  )
}

checks <- rbind(

  check_num(
    "Analytic n",
    expected$n,
    analytic_n,
    0
  ),

  check_num(
    "Survey df",
    expected$design_df,
    design_df,
    0
  ),

  check_num(
    "Interaction F",
    expected$interaction_F,
    interaction_F
  ),

  check_num(
    "Interaction p",
    expected$interaction_p,
    interaction_p
  ),

  check_num(
    "Female slope",
    expected$female_beta,
    female$beta
  ),

  check_num(
    "Female SE",
    expected$female_se,
    female$SE
  ),

  check_num(
    "Female CI low",
    expected$female_ci_low,
    female$CI_low
  ),

  check_num(
    "Female CI high",
    expected$female_ci_high,
    female$CI_high
  ),

  check_num(
    "Female p",
    expected$female_p,
    female$p
  ),

  check_num(
    "Male slope",
    expected$male_beta,
    male$beta
  ),

  check_num(
    "Male SE",
    expected$male_se,
    male$SE
  ),

  check_num(
    "Male CI low",
    expected$male_ci_low,
    male$CI_low
  ),

  check_num(
    "Male CI high",
    expected$male_ci_high,
    male$CI_high
  ),

  check_num(
    "Male p",
    expected$male_p,
    male$p
  ),

  check_num(
    "Male-minus-female contrast",
    expected$male_minus_female,
    difference$beta
  ),

  check_num(
    "Difference CI low",
    expected$difference_ci_low,
    difference$CI_low
  ),

  check_num(
    "Difference CI high",
    expected$difference_ci_high,
    difference$CI_high
  )
)

cat("PYTHON-vs-R LOCKED VALIDATION\n")

print(
  checks,
  digits = 12,
  row.names = FALSE
)

cat("\n")

# ------------------------------------------------------------
# Final validation status
# ------------------------------------------------------------

if (all(checks$status == "PASS")) {

  cat("============================================================\n")
  cat("STUDY #2 PRIMARY CROSS-SOFTWARE VALIDATION: PASS\n")
  cat("============================================================\n")
  cat("\n")
  cat("R independently reproduced:\n")
  cat("  - the locked analytic cohort\n")
  cat("  - survey design degrees of freedom\n")
  cat("  - the HOMA-IR x sex interaction\n")
  cat("  - the female HOMA-IR slope\n")
  cat("  - the male HOMA-IR slope\n")
  cat("  - the male-minus-female slope contrast and 95% CI\n")
  cat("\n")

  quit(
    status = 0
  )

} else {

  cat("============================================================\n")
  cat("STUDY #2 CROSS-SOFTWARE VALIDATION: DISCREPANCY DETECTED\n")
  cat("============================================================\n")
  cat("\n")
  cat("Do not tune the R analysis to force agreement.\n")
  cat("Inspect the failed target(s) before changing any code.\n")
  cat("\n")

  quit(
    status = 2
  )
}
