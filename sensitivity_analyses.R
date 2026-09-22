# Sensitivity analyses for primary outcome (SATE)
## All persons with missing data have only failure or only success out of 15 answers

panel_test <- readRDS("panel_test.rds")

# All with missing data is coded with failure on every answer
sens_0 <- panel_test |>
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

range(sens_0$scenario_sum)
range(sens_0$scenario_fail)

sens_0 |>
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

sens0_rr <- glm(
  cbind(scenario_sum, scenario_fail) ~ arm,
  family = binomial(link = "log"),
  data = sens_0
)

summary(sens0_rr)

# number of correct answers in each arm
sens_0 |>
  group_by(arm) |>
  summarise(
    correct = sum(scenario_sum),
    not_corr = sum(scenario_fail),
    total = sum(correct + not_corr),
    pct = 100 * correct / total
  )

# -------------------------------------------------------------------------
# 3. Check for over-dispersion (Pearson chi-square / residual df)
# -------------------------------------------------------------------------

dispersion <- sum(residuals(sens0_rr, type = "pearson")^2) / df.residual(sens0_rr)
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
  V <- vcovHC(sens0_rr, type = "HC3")
} else {
  message(sprintf(
    "Dispersion = %.2f (<= 1.2): no evidence of over-dispersion, using model-based SEs.",
    dispersion
  ))
  V <- vcov(sens0_rr)
}

coeftest(sens0_rr, vcov = V)

# -------------------------------------------------------------------------
# 4. Bonferroni adjustment for the 3 pairwise comparisons vs. V1_control
# -------------------------------------------------------------------------

beta <- coef(sens0_rr)[-1]          # drop intercept
se   <- sqrt(diag(V))[-1]

m <- length(beta)                 # number of comparisons vs. reference (3)
alpha_bonf <- 0.05 / m
z_crit <- qnorm(1 - alpha_bonf / 2)

z_value <- beta / se
p_raw   <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
p_bonf  <- p.adjust(p_raw, method = "bonferroni")

## Risk difference (RD), re-expressed from the RR.
## RD is a non-linear function of two correlated model parameters -- the
## control-arm log-risk (intercept, b0) and the log-RR (b1) -- so a delta
## method is used to propagate uncertainty in BOTH into the SE of RD, rather
## than treating the control-arm risk as fixed:
##   p0 = exp(b0), p1 = exp(b0 + b1), RD = p1 - p0
##   d(RD)/d(b0) = p1 - p0 = RD ; d(RD)/d(b1) = p1
##   Var(RD) = (dRD/db0)^2 Var(b0) + (dRD/db1)^2 Var(b1) + 2 (dRD/db0)(dRD/db1) Cov(b0,b1)
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

coefs_rr <- coef(sens0_rr)
rd_sate <- purrr::map_dfr(seq_along(beta) + 1, function(i) {
  out <- delta_rd(coefs_rr, V, i, ref_idx = 1)
  tibble(term = names(coefs_rr)[i], RD = out["rd"], RD_SE = out["se"])
}) |>
  mutate(RD_CI_lower = RD - z_crit * RD_SE, RD_CI_upper = RD + z_crit * RD_SE)

resultat <- tibble(
  term      = names(beta),
  RR        = exp(beta),
  SE        = if (use_robust_se) "robust (HC3)" else "model-based",
  CI_lower  = exp(beta - z_crit * se),
  CI_upper  = exp(beta + z_crit * se),
  p_bonferroni = p_bonf
) |>
  left_join(rd_sate, by = "term") |>
  mutate(
    Comparison = recode(
      term,
      "armV2_sentence"             = "V2 (sentence) vs. V1 (control)",
      "armV3_definitions"          = "V3 (definitions) vs. V1 (control)",
      "armV4_sentence_definitions" = "V4 (sentence + definitions) vs. V1 (control)"
    ),
    `RR (95% Bonferroni-adjusted CI)` = sprintf("%.2f (%.2f to %.2f)", RR, CI_lower, CI_upper),
    `RD (95% Bonferroni-adjusted CI, pct. points)` = sprintf(
      "%.1f (%.1f to %.1f)", 100 * RD, 100 * RD_CI_lower, 100 * RD_CI_upper
    ),
    `Bonferroni-adjusted p value` = case_when(
      p_bonferroni < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_bonferroni)
    )
  ) |>
  select(Comparison, `RR (95% Bonferroni-adjusted CI)`,
         `RD (95% Bonferroni-adjusted CI, pct. points)`, `Bonferroni-adjusted p value`)

resultat

#-------------------------------------------------------------------------------

panel_test <- readRDS("panel_test.rds")

# All with missing data is coded with failure on every answer
sens_1 <- panel_test |>
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
  
  # Change missing to only success
  mutate(sc1_1 = if_else(is.na(answer_time_ms), 1, sc1_1),
         sc1_2 = if_else(is.na(answer_time_ms), 1, sc1_2),
         sc1_3 = if_else(is.na(answer_time_ms), 1, sc1_3),
         
         sc2_1 = if_else(is.na(answer_time_ms), 1, sc2_1),
         sc2_2 = if_else(is.na(answer_time_ms), 1, sc2_2),
         sc2_3 = if_else(is.na(answer_time_ms), 1, sc2_3),
         
         sc3_1 = if_else(is.na(answer_time_ms), 1, sc3_1),
         sc3_2 = if_else(is.na(answer_time_ms), 1, sc3_2),
         sc3_3 = if_else(is.na(answer_time_ms), 1, sc3_3),
         
         sc4_1 = if_else(is.na(answer_time_ms), 1, sc4_1),
         sc4_2 = if_else(is.na(answer_time_ms), 1, sc4_2),
         sc4_3 = if_else(is.na(answer_time_ms), 1, sc4_3),
         
         sc5_1 = if_else(is.na(answer_time_ms), 1, sc5_1),
         sc5_2 = if_else(is.na(answer_time_ms), 1, sc5_2),
         sc5_3 = if_else(is.na(answer_time_ms), 1, sc5_3)) |>
  
  # Sum of all 15 items
  mutate(scenario_sum = rowSums(across(c(sc1_1, sc1_2, sc1_3,
                                         sc2_1, sc2_2, sc2_3,
                                         sc3_1, sc3_2, sc3_3,
                                         sc4_1, sc4_2, sc4_3,
                                         sc5_1, sc5_2, sc5_3)))) |>
  mutate(scenario_fail = 15 - scenario_sum) |>
  # Reference arm for all comparisons
  mutate(arm = relevel(factor(arm), ref = "V1_control"))

range(sens_1$scenario_sum)
range(sens_1$scenario_fail)

sens_1 |>
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

sens1_rr <- glm(
  cbind(scenario_sum, scenario_fail) ~ arm,
  family = binomial(link = "log"),
  data = sens_1
)

summary(sens1_rr)

# number of correct answers in each arm
sens_1 |>
  group_by(arm) |>
  summarise(
    correct = sum(scenario_sum),
    not_corr = sum(scenario_fail),
    total = sum(correct + not_corr),
    pct = 100 * correct / total
  )

# -------------------------------------------------------------------------
# 3. Check for over-dispersion (Pearson chi-square / residual df)
# -------------------------------------------------------------------------

dispersion <- sum(residuals(sens1_rr, type = "pearson")^2) / df.residual(sens1_rr)
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
  V <- vcovHC(sens1_rr, type = "HC3")
} else {
  message(sprintf(
    "Dispersion = %.2f (<= 1.2): no evidence of over-dispersion, using model-based SEs.",
    dispersion
  ))
  V <- vcov(sens1_rr)
}

coeftest(sens1_rr, vcov = V)

# -------------------------------------------------------------------------
# 4. Bonferroni adjustment for the 3 pairwise comparisons vs. V1_control
# -------------------------------------------------------------------------

beta <- coef(sens1_rr)[-1]          # drop intercept
se   <- sqrt(diag(V))[-1]

m <- length(beta)                 # number of comparisons vs. reference (3)
alpha_bonf <- 0.05 / m
z_crit <- qnorm(1 - alpha_bonf / 2)

z_value <- beta / se
p_raw   <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
p_bonf  <- p.adjust(p_raw, method = "bonferroni")

## Risk difference (RD), re-expressed from the RR.
## RD is a non-linear function of two correlated model parameters -- the
## control-arm log-risk (intercept, b0) and the log-RR (b1) -- so a delta
## method is used to propagate uncertainty in BOTH into the SE of RD, rather
## than treating the control-arm risk as fixed:
##   p0 = exp(b0), p1 = exp(b0 + b1), RD = p1 - p0
##   d(RD)/d(b0) = p1 - p0 = RD ; d(RD)/d(b1) = p1
##   Var(RD) = (dRD/db0)^2 Var(b0) + (dRD/db1)^2 Var(b1) + 2 (dRD/db0)(dRD/db1) Cov(b0,b1)
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

coefs_rr <- coef(sens1_rr)
rd_sate <- purrr::map_dfr(seq_along(beta) + 1, function(i) {
  out <- delta_rd(coefs_rr, V, i, ref_idx = 1)
  tibble(term = names(coefs_rr)[i], RD = out["rd"], RD_SE = out["se"])
}) |>
  mutate(RD_CI_lower = RD - z_crit * RD_SE, RD_CI_upper = RD + z_crit * RD_SE)

resultat <- tibble(
  term      = names(beta),
  RR        = exp(beta),
  SE        = if (use_robust_se) "robust (HC3)" else "model-based",
  CI_lower  = exp(beta - z_crit * se),
  CI_upper  = exp(beta + z_crit * se),
  p_bonferroni = p_bonf
) |>
  left_join(rd_sate, by = "term") |>
  mutate(
    Comparison = recode(
      term,
      "armV2_sentence"             = "V2 (sentence) vs. V1 (control)",
      "armV3_definitions"          = "V3 (definitions) vs. V1 (control)",
      "armV4_sentence_definitions" = "V4 (sentence + definitions) vs. V1 (control)"
    ),
    `RR (95% Bonferroni-adjusted CI)` = sprintf("%.2f (%.2f to %.2f)", RR, CI_lower, CI_upper),
    `RD (95% Bonferroni-adjusted CI, pct. points)` = sprintf(
      "%.1f (%.1f to %.1f)", 100 * RD, 100 * RD_CI_lower, 100 * RD_CI_upper
    ),
    `Bonferroni-adjusted p value` = case_when(
      p_bonferroni < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_bonferroni)
    )
  ) |>
  select(Comparison, `RR (95% Bonferroni-adjusted CI)`,
         `RD (95% Bonferroni-adjusted CI, pct. points)`, `Bonferroni-adjusted p value`)

resultat

#-------------------------------------------------------------------------------

