# Subgroup analyses for primary outcome
## One analysis performed as described in prim_outcome.R (log-binomial RR vs.
## V1_control, Bonferroni-adjusted for 3 comparisons), for each level of the
## following subgroups. Mapping decisions agreed with study team (native
## panel_test categories don't align exactly with the originally-planned
## bins -- see notes at each variable below):
#
# Age group (native age_group has 6x 10-year bins from 16; collapsed to 3):
##    <35    (16-24, 25-34)
##    35-54  (35-44, 45-54)
##    55+    (55-64, 65 or older)
# Gender (native gender also has Other/Prefer not to say, n = 2/4 -- too few
## for a separate stratum, so this subgroup analysis is restricted to):
##    Women
##    Men
# Education level (native education has 5 levels; collapsed to 3):
##    Low    (Primary school or lower)
##    Middle (Upper secondary school, Vocational college)
##    High   (University/college, <=4 years or >4 years)

# ------------------------------------------------------------------------------
# Start by running descreptive_analysis.R (to start a new result folder) and
## run prim_outcome.R (to ensure correct variables in the data frame panel_test)

# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(sandwich)
library(lmtest)
library(gt)

required_cols <- c("scenario_sum", "scenario_fail", "arm", "age_group", "gender", "education")
missing_cols <- setdiff(required_cols, names(panel_test))
if (length(missing_cols) > 0) {
  stop(
    "panel_test is missing: ", paste(missing_cols, collapse = ", "),
    ". Run prim_outcome.R (section 1) first to build the outcome variables."
  )
}

# -------------------------------------------------------------------------
# 1. Define subgroup variables
# -------------------------------------------------------------------------

panel_sub <- panel_test |>
  mutate(
    age_group_3 = case_when(
      age_group %in% c("16-24", "25-34")        ~ "<35",
      age_group %in% c("35-44", "45-54")        ~ "35-54",
      age_group %in% c("55-64", "65 or older")  ~ "55+",
      TRUE ~ NA_character_
    ) |> factor(levels = c("<35", "35-54", "55+")),

    gender_2 = case_when(
      gender == "Female" ~ "Women",
      gender == "Male"   ~ "Men",
      TRUE ~ NA_character_          # Other / Prefer not to say -> excluded
    ) |> factor(levels = c("Women", "Men")),

    education_3 = case_when(
      education == "Primary school or lower" ~ "Low",
      education %in% c("Upper secondary school", "Vocational college") ~ "Middle",
      education %in% c("University/college, 4 years or less",
                        "University/college, more than 4 years") ~ "High",
      TRUE ~ NA_character_
    ) |> factor(levels = c("Low", "Middle", "High"))
  )

# Cell counts per subgroup level x arm -- inspect for strata too sparse for a
# stable log-binomial fit before trusting the results below.
panel_sub |> count(age_group_3, arm) |> pivot_wider(names_from = arm, values_from = n, values_fill = 0)
panel_sub |> count(gender_2, arm)    |> pivot_wider(names_from = arm, values_from = n, values_fill = 0)
panel_sub |> count(education_3, arm) |> pivot_wider(names_from = arm, values_from = n, values_fill = 0)

n_gender_excluded <- sum(is.na(panel_sub$gender_2))
if (n_gender_excluded > 0) {
  message(sprintf(
    "%d respondent(s) with gender = Other/Prefer not to say excluded from the gender subgroup analysis.",
    n_gender_excluded
  ))
}

# -------------------------------------------------------------------------
# 2. Reusable RR/RD analysis, mirroring prim_outcome.R sections 2-4, applied
#    within a single subgroup level (one arm-comparison model per level)
# -------------------------------------------------------------------------

## Delta-method SE for the risk difference (RD), re-expressed from the RR --
## see prim_outcome.R section 4 for the full derivation/rationale.
delta_rd <- function(coefs, Vmat, term_idx, ref_idx = 1) {
  b0 <- unname(coefs[ref_idx]); b1 <- unname(coefs[term_idx])
  v00 <- Vmat[ref_idx, ref_idx]
  v11 <- Vmat[term_idx, term_idx]
  v01 <- Vmat[ref_idx, term_idx]
  p0 <- exp(b0)
  p1 <- exp(b0 + b1)
  rd <- p1 - p0
  d_b0 <- rd
  d_b1 <- p1
  var_rd <- d_b0^2 * v00 + d_b1^2 * v11 + 2 * d_b0 * d_b1 * v01
  c(rd = unname(rd), se = unname(sqrt(var_rd)))
}

arm_labels_recode <- c(
  "armV2_sentence"             = "V2 (sentence) vs. V1 (control)",
  "armV3_definitions"          = "V3 (definitions) vs. V1 (control)",
  "armV4_sentence_definitions" = "V4 (sentence + definitions) vs. V1 (control)"
)

run_subgroup_rr <- function(data, subgroup_label) {
  data <- data |> mutate(arm = droplevels(arm))

  if (nlevels(data$arm) < 2 || any(table(data$arm) == 0)) {
    warning(sprintf(
      "Skipping subgroup '%s': fewer than 2 arms with data (n = %d).",
      subgroup_label, nrow(data)
    ))
    return(NULL)
  }

  fit <- glm(
    cbind(scenario_sum, scenario_fail) ~ arm,
    family = binomial(link = "log"),
    data = data
  )

  dispersion <- sum(residuals(fit, type = "pearson")^2) / df.residual(fit)
  use_robust_se <- dispersion > 1.2
  V <- if (use_robust_se) vcovHC(fit, type = "HC3") else vcov(fit)

  beta <- coef(fit)[-1]
  se   <- sqrt(diag(V))[-1]

  m <- length(beta)
  alpha_bonf <- 0.05 / m
  z_crit <- qnorm(1 - alpha_bonf / 2)

  z_value <- beta / se
  p_raw   <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
  p_bonf  <- p.adjust(p_raw, method = "bonferroni")

  coefs <- coef(fit)
  rd <- purrr::map_dfr(seq_along(beta) + 1, function(i) {
    out <- delta_rd(coefs, V, i, ref_idx = 1)
    tibble(term = names(coefs)[i], RD = out["rd"], RD_SE = out["se"])
  }) |>
    mutate(RD_CI_lower = RD - z_crit * RD_SE, RD_CI_upper = RD + z_crit * RD_SE)

  tibble(
    term         = names(beta),
    n            = nrow(data),
    RR           = exp(beta),
    CI_lower     = exp(beta - z_crit * se),
    CI_upper     = exp(beta + z_crit * se),
    p_bonferroni = p_bonf,
    dispersion   = dispersion,
    SE_type      = if (use_robust_se) "robust (HC3)" else "model-based"
  ) |>
    left_join(rd, by = "term") |>
    mutate(
      Subgroup   = subgroup_label,
      Comparison = recode(term, !!!arm_labels_recode),
      `RR (95% Bonferroni-adjusted CI)` = sprintf("%.2f (%.2f to %.2f)", RR, CI_lower, CI_upper),
      `RD (95% Bonferroni-adjusted CI, pct. points)` = sprintf(
        "%.1f (%.1f to %.1f)", 100 * RD, 100 * RD_CI_lower, 100 * RD_CI_upper
      ),
      `Bonferroni-adjusted p value` = case_when(
        p_bonferroni < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p_bonferroni)
      )
    ) |>
    select(Subgroup, Comparison, n, `RR (95% Bonferroni-adjusted CI)`,
           `RD (95% Bonferroni-adjusted CI, pct. points)`, `Bonferroni-adjusted p value`,
           dispersion, SE_type)
}

## Runs run_subgroup_rr() separately within each observed level of group_var
## (NAs in group_var, e.g. excluded gender categories, are dropped first).
run_subgroup_dimension <- function(data, group_var) {
  data <- data |> filter(!is.na(.data[[group_var]]))
  splits <- split(data, data[[group_var]])
  splits <- splits[sapply(splits, nrow) > 0]
  map_dfr(names(splits), function(lvl) run_subgroup_rr(splits[[lvl]], lvl))
}

# -------------------------------------------------------------------------
# 3. Run each subgroup dimension
# -------------------------------------------------------------------------

resultat_age       <- run_subgroup_dimension(panel_sub, "age_group_3")
resultat_gender     <- run_subgroup_dimension(panel_sub, "gender_2")
resultat_education <- run_subgroup_dimension(panel_sub, "education_3")

resultat_age
resultat_gender
resultat_education

resultat_subgroups <- bind_rows(
  mutate(resultat_age, Dimension = "Age group", .before = Subgroup),
  mutate(resultat_gender, Dimension = "Gender", .before = Subgroup),
  mutate(resultat_education, Dimension = "Education level", .before = Subgroup)
)

resultat_subgroups

write.csv(
  resultat_subgroups |> select(-dispersion, -SE_type),
  file.path(results_dir, "table6_subgroup_analyses.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 4. Publication-ready table (gt)
# -------------------------------------------------------------------------

subgroup_gt <- resultat_subgroups |>
  select(Dimension, Subgroup, Comparison, n,
         `RR (95% Bonferroni-adjusted CI)`,
         `RD (95% Bonferroni-adjusted CI, pct. points)`,
         `Bonferroni-adjusted p value`) |>
  gt(groupname_col = "Dimension", rowname_col = "Subgroup") |>
  tab_header(
    title = "Table 6. Subgroup analyses of the primary outcome",
    subtitle = "Log-binomial regression vs. V1 (control), within each subgroup level; Bonferroni-adjusted for 3 comparisons per level"
  ) |>
  cols_label(
    Comparison = "Comparison",
    n = "n",
    `RR (95% Bonferroni-adjusted CI)` = "RR (95% Bonferroni-adjusted CI)",
    `RD (95% Bonferroni-adjusted CI, pct. points)` = "RD, pct. points (95% Bonferroni-adjusted CI)",
    `Bonferroni-adjusted p value` = "Bonferroni-adjusted p value"
  ) |>
  cols_align(align = "left", columns = Comparison) |>
  cols_align(align = "center", columns = c(n, `RR (95% Bonferroni-adjusted CI)`,
                                            `RD (95% Bonferroni-adjusted CI, pct. points)`,
                                            `Bonferroni-adjusted p value`)) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) |>
  tab_footnote(
    footnote = paste(
      "Age group: native 6-level age_group collapsed to <35 (16-24, 25-34),",
      "35-54 (35-44, 45-54), 55+ (55-64, 65 or older). Gender: restricted to",
      "Female/Male (Other/Prefer not to say excluded, n too small for a",
      "separate stratum). Education: native 5-level education collapsed to",
      "Low (primary school or lower), Middle (upper secondary school,",
      "vocational college), High (university/college, any duration)."
    )
  ) |>
  tab_footnote(
    footnote = paste(
      "Each row is a separate model fit within that subgroup level only",
      "(not adjusted for other subgroups or interaction-tested against the",
      "overall trial estimate); small strata (e.g. age <35, education Low)",
      "have wide confidence intervals and should be interpreted cautiously.",
      "SE type (model-based binomial vs. robust HC3, chosen per stratum by",
      "the same dispersion > 1.2 rule as prim_outcome.R) is not shown but",
      "is available in resultat_subgroups$SE_type."
    )
  ) |>
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    column_labels.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.bottom.style = "solid"
  ) |>
  opt_table_font(font = "Times New Roman")

subgroup_gt

gtsave_safe <- function(data, path) {
  tmp <- tempfile(fileext = paste0(".", tools::file_ext(path)))
  gtsave(data, tmp)
  file.copy(tmp, path, overwrite = TRUE)
  file.remove(tmp)
  invisible(path)
}

gtsave_safe(subgroup_gt, file.path(results_dir, "table6_subgroup_analyses.docx"))
gtsave_safe(subgroup_gt, file.path(results_dir, "table6_subgroup_analyses.rtf"))
gtsave_safe(subgroup_gt, file.path(results_dir, "table6_subgroup_analyses.html"))

#-------------------------------------------------------------------------------
