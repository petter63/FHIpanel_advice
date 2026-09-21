# Primary outcome
## Number of correct answers across the 5 scenarios x 3 items (0-15),
## analysed as a relative risk (RR) via log-binomial regression.

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(sandwich)
library(lmtest)
library(survey)
library(gt)

# -------------------------------------------------------------------------
# 1. Code correct answers and build the outcome (successes / failures out of 15)
# -------------------------------------------------------------------------

panel_test <- panel_test |>
  mutate(sc1_1 = if_else(scenario1_item1 %in% c("Likely", "Very likely"), 1, 0),
         sc1_2 = if_else(scenario1_item2_careful %in% c("Likely", "Very likely"), 1, 0),
         sc1_3 = if_else(scenario1_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0),

         sc2_1 = if_else(scenario2_item1 %in% c("Likely", "Very likely"), 1, 0),
         sc2_2 = if_else(scenario2_item2_careful %in% c("Likely", "Very likely"), 1, 0),
         sc2_3 = if_else(scenario2_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0),

         sc3_1 = if_else(scenario3_item1 %in% c("Very unlikely", "Unlikely"), 1, 0),
         sc3_2 = if_else(scenario3_item2_careful %in% c("Likely", "Very likely"), 1, 0),
         sc3_3 = if_else(scenario3_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0),

         sc4_1 = if_else(scenario4_item1 %in% c("Very unlikely", "Unlikely"), 1, 0),
         sc4_2 = if_else(scenario4_item2_careful %in% c("Likely", "Very likely"), 1, 0),
         sc4_3 = if_else(scenario4_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0),

         sc5_1 = if_else(scenario5_item1 %in% c("Very unlikely", "Unlikely"), 1, 0),
         sc5_2 = if_else(scenario5_item2_careful %in% c("Likely", "Very likely"), 1, 0),
         sc5_3 = if_else(scenario5_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0)) |>

  # Sum of all 15 items
  mutate(scenario_sum = rowSums(across(c(sc1_1, sc1_2, sc1_3,
                                          sc2_1, sc2_2, sc2_3,
                                          sc3_1, sc3_2, sc3_3,
                                          sc4_1, sc4_2, sc4_3,
                                          sc5_1, sc5_2, sc5_3)))) |>
  mutate(scenario_fail = 15 - scenario_sum) |>
  # Reference arm for all comparisons
  mutate(arm = relevel(factor(arm), ref = "V1_control"))

range(panel_test$scenario_sum)
range(panel_test$scenario_fail)

panel_test |>
  group_by(arm) |>
  summarise(
    mean   = mean(scenario_sum, na.rm = TRUE),
    sd     = sd(scenario_sum, na.rm = TRUE),
    median = median(scenario_sum, na.rm = TRUE),
    IQR    = IQR(scenario_sum, na.rm = TRUE)
  )

# -------------------------------------------------------------------------
# 2. Log-binomial regression (RR of a correct answer, vs. V1_control)
# -------------------------------------------------------------------------

fit_rr <- glm(
  cbind(scenario_sum, scenario_fail) ~ arm,
  family = binomial(link = "log"),
  data = panel_test
)

summary(fit_rr)

# -------------------------------------------------------------------------
# 3. Check for over-dispersion (Pearson chi-square / residual df)
# -------------------------------------------------------------------------

dispersion <- sum(residuals(fit_rr, type = "pearson")^2) / df.residual(fit_rr)
dispersion

# Rule of thumb: dispersion appreciably above 1 indicates the binomial
# variance assumption doesn't hold, so we fall back to a robust
# (heteroscedasticity-consistent) sandwich estimator for the SEs.
use_robust_se <- dispersion > 1.2

if (use_robust_se) {
  message(sprintf(
    "Dispersion = %.2f (> 1.2): using HC3 robust SEs.",
    dispersion
  ))
  V <- vcovHC(fit_rr, type = "HC3")
} else {
  message(sprintf(
    "Dispersion = %.2f (<= 1.2): no evidence of over-dispersion, using model-based SEs.",
    dispersion
  ))
  V <- vcov(fit_rr)
}

coeftest(fit_rr, vcov = V)

# -------------------------------------------------------------------------
# 4. Bonferroni adjustment for the 3 pairwise comparisons vs. V1_control
# -------------------------------------------------------------------------

beta <- coef(fit_rr)[-1]          # drop intercept
se   <- sqrt(diag(V))[-1]

m <- length(beta)                 # number of comparisons vs. reference (3)
alpha_bonf <- 0.05 / m
z_crit <- qnorm(1 - alpha_bonf / 2)

z_value <- beta / se
p_raw   <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
p_bonf  <- p.adjust(p_raw, method = "bonferroni")

resultat <- tibble(
  term      = names(beta),
  RR        = exp(beta),
  SE        = if (use_robust_se) "robust (HC3)" else "model-based",
  CI_lower  = exp(beta - z_crit * se),
  CI_upper  = exp(beta + z_crit * se),
  p_bonferroni = p_bonf
) |>
  mutate(
    Comparison = recode(
      term,
      "armV2_sentence"             = "V2 (sentence) vs. V1 (control)",
      "armV3_definitions"          = "V3 (definitions) vs. V1 (control)",
      "armV4_sentence_definitions" = "V4 (sentence + definitions) vs. V1 (control)"
    ),
    `RR (95% Bonferroni-adjusted CI)` = sprintf("%.2f (%.2f to %.2f)", RR, CI_lower, CI_upper),
    `Bonferroni-adjusted p value` = case_when(
      p_bonferroni < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_bonferroni)
    )
  ) |>
  select(Comparison, `RR (95% Bonferroni-adjusted CI)`, `Bonferroni-adjusted p value`)

resultat

# -------------------------------------------------------------------------
# 5. Publication-ready table (gt)
# -------------------------------------------------------------------------

resultat_gt <- resultat |>
  gt() |>
  tab_header(
    title = "Table 3. Relative risk of a correct comprehension answer",
    subtitle = sprintf(
      "Log-binomial regression vs. V1 (control); %s; Bonferroni-adjusted for %d comparisons",
      if (use_robust_se) "robust (HC3) SEs used due to over-dispersion" else "model-based SEs",
      m
    )
  ) |>
  cols_label(
    Comparison = "Comparison",
    `RR (95% Bonferroni-adjusted CI)` = "RR (95% Bonferroni-adjusted CI)",
    `Bonferroni-adjusted p value` = "Bonferroni-adjusted p value"
  ) |>
  cols_align(align = "left", columns = Comparison) |>
  cols_align(align = "center", columns = c(`RR (95% Bonferroni-adjusted CI)`, `Bonferroni-adjusted p value`)) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_footnote(
    footnote = sprintf(
      "Outcome: number of correct answers out of 15 comprehension items. Dispersion statistic = %.2f.",
      dispersion
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

resultat_gt

# Saving straight into a OneDrive-synced folder occasionally makes pandoc
# fail with "error 22" because OneDrive briefly locks the file mid-write.
# Writing to a local temp file first and copying it into place avoids this.
gtsave_safe <- function(data, path) {
  tmp <- tempfile(fileext = paste0(".", tools::file_ext(path)))
  gtsave(data, tmp)
  file.copy(tmp, path, overwrite = TRUE)
  file.remove(tmp)
  invisible(path)
}

gtsave_safe(resultat_gt, file.path(results_dir, "table3_primary_outcome_RR.docx"))
gtsave_safe(resultat_gt, file.path(results_dir, "table3_primary_outcome_RR.rtf"))
gtsave_safe(resultat_gt, file.path(results_dir, "table3_primary_outcome_RR.html"))

#-------------------------------------------------------------------------------
# 6. Population average treatment effect (PATE), using post-stratification
#    weights from weighting.R (panel_weighting.rds)
# -------------------------------------------------------------------------
## The RR above is a SAMPLE average treatment effect: it describes the
## effect for respondents as recruited, whose age/gender/education mix may
## not match the general population. Re-weighting each arm to the national
## population (age, gender, education; see weighting.R) targets the
## population average treatment effect (PATE) instead -- i.e. the effect we
## would expect if the trial's population had matched Norway's adult
## population on these variables.
##
## Modelling choices, mirroring section 2 as closely as the weighting
## allows:
##  - Item-level (long) format: one row per respondent x comprehension item
##    (15 rows/respondent) with a 0/1 `correct` outcome, rather than
##    aggregated cbind(scenario_sum, scenario_fail) counts -- this keeps the
##    post-stratification weight attached at the respondent level (it would
##    otherwise be ambiguous whether a weight applies to the count or the
##    denominator) and lets `ids = ~participant_id` account for the
##    within-respondent correlation across the 15 items.
##  - `svyglm()` (not `glm()`) so that the standard errors are the
##    linearised, design-based ones appropriate for weighted/clustered data,
##    rather than model-based binomial SEs.
##  - quasibinomial(link = "log") for the same RR interpretation as
##    section 2's log-binomial model.
##  - Respondents excluded from raking (missing age/gender/education, or
##    gender = Other/Prefer not to say -- see weighting.R) keep weight = 1,
##    i.e. contribute to the PATE estimate unadjusted, since there is no
##    population benchmark to re-weight them against.

panel_weighting <- readRDS("panel_weighting.rds") |>
  select(participant_id, weight)

panel_test_w <- panel_test |>
  left_join(panel_weighting, by = "participant_id")

n_unmatched <- sum(is.na(panel_test_w$weight))
if (n_unmatched > 0) {
  warning(n_unmatched, " respondent(s) in panel_test have no matching weight ",
          "in panel_weighting.rds -- check that weighting.R was run on the ",
          "same panel_test.rds.")
}

itemized_w <- panel_test_w |>
  select(participant_id, arm, weight,
         sc1_1, sc1_2, sc1_3, sc2_1, sc2_2, sc2_3, sc3_1, sc3_2, sc3_3,
         sc4_1, sc4_2, sc4_3, sc5_1, sc5_2, sc5_3) |>
  pivot_longer(cols = starts_with("sc"), names_to = "item", values_to = "correct")

design_w <- svydesign(ids = ~participant_id, weights = ~weight, data = itemized_w)

fit_pate <- svyglm(correct ~ arm, design = design_w, family = quasibinomial(link = "log"))

summary(fit_pate)

# -------------------------------------------------------------------------
# 6b. Bonferroni adjustment for the 3 pairwise comparisons vs. V1_control
# -------------------------------------------------------------------------

beta_pate <- coef(fit_pate)[-1]
se_pate   <- sqrt(diag(vcov(fit_pate)))[-1]  # already design-based/robust via svyglm

m_pate      <- length(beta_pate)
alpha_bonf_pate <- 0.05 / m_pate
# svyglm uses a t-reference distribution (df = number of clusters - 1)
df_pate     <- df.residual(fit_pate)
t_crit_pate <- qt(1 - alpha_bonf_pate / 2, df = df_pate)

t_value_pate <- beta_pate / se_pate
p_raw_pate   <- 2 * pt(abs(t_value_pate), df = df_pate, lower.tail = FALSE)
p_bonf_pate  <- p.adjust(p_raw_pate, method = "bonferroni")

resultat_pate <- tibble(
  term      = names(beta_pate),
  RR        = exp(beta_pate),
  CI_lower  = exp(beta_pate - t_crit_pate * se_pate),
  CI_upper  = exp(beta_pate + t_crit_pate * se_pate),
  p_bonferroni = p_bonf_pate
) |>
  mutate(
    Comparison = recode(
      term,
      "armV2_sentence"             = "V2 (sentence) vs. V1 (control)",
      "armV3_definitions"          = "V3 (definitions) vs. V1 (control)",
      "armV4_sentence_definitions" = "V4 (sentence + definitions) vs. V1 (control)"
    ),
    `PATE RR (95% Bonferroni-adjusted CI)` = sprintf("%.2f (%.2f to %.2f)", RR, CI_lower, CI_upper),
    `Bonferroni-adjusted p value` = case_when(
      p_bonferroni < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_bonferroni)
    )
  ) |>
  select(Comparison, `PATE RR (95% Bonferroni-adjusted CI)`, `Bonferroni-adjusted p value`)

resultat_pate

# -------------------------------------------------------------------------
# 6c. Publication-ready table (gt)
# -------------------------------------------------------------------------

resultat_pate_gt <- resultat_pate |>
  gt() |>
  tab_header(
    title = "Table 4. Population average treatment effect (PATE) on a correct comprehension answer",
    subtitle = sprintf(
      paste(
        "Item-level log-binomial (quasibinomial) regression vs. V1 (control),",
        "post-stratification weighted (age, gender, education) and clustered",
        "by participant; design-based SEs; Bonferroni-adjusted for %d comparisons"
      ),
      m_pate
    )
  ) |>
  cols_label(
    Comparison = "Comparison",
    `PATE RR (95% Bonferroni-adjusted CI)` = "PATE RR (95% Bonferroni-adjusted CI)",
    `Bonferroni-adjusted p value` = "Bonferroni-adjusted p value"
  ) |>
  cols_align(align = "left", columns = Comparison) |>
  cols_align(align = "center", columns = c(`PATE RR (95% Bonferroni-adjusted CI)`, `Bonferroni-adjusted p value`)) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_footnote(
    footnote = paste(
      "Weights from panel_weighting.rds (weighting.R): raked separately within",
      "each arm to national population margins for age group, gender, and",
      "education level (SSB). Respondents with no population benchmark",
      "(missing demographics, or gender = Other/Prefer not to say) retain",
      "weight = 1."
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

resultat_pate_gt

gtsave_safe(resultat_pate_gt, file.path(results_dir, "table4_primary_outcome_PATE.docx"))
gtsave_safe(resultat_pate_gt, file.path(results_dir, "table4_primary_outcome_PATE.rtf"))
gtsave_safe(resultat_pate_gt, file.path(results_dir, "table4_primary_outcome_PATE.html"))

#-------------------------------------------------------------------------------

