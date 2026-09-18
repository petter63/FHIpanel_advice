# Primary outcome
## Number of correct answers across the 5 scenarios x 3 items (0-15),
## analysed as a relative risk (RR) via log-binomial regression.

library(dplyr)
library(purrr)
library(broom)
library(sandwich)
library(lmtest)
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
