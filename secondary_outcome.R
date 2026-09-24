# Secondary analyses
## Final count: 5 research questions with a dedicated analysis (SQ1, SQ2,
## SQ3, SQ5; SQ4 has no separate analysis -- it is interpreted from the
## SQ1-SQ3 results, per the protocol note below) and 12 secondary analyses
## in total: SO1_1, SO1_2-SO1_6 (6), SO2_1 (1), SO3_1-SO3_2 (2),
## SO5_1-SO5_3 (3).

# --------------------------------------------------------------------
# SQ1. Does providing additional explanations of key terms used in infection
## control advice increase the likelihood that members of the public report 
## intended adherence behaviour with symptoms of respiratory infection?

# SO1_1. Meeting or exceeding threshold understanding for intended behaviours.
## Outcome: Responding correctly for >= 12 (80%) of the 15 questions for the 5
## hypothetical scenarios
## Dichotomt outcome for each person: scenario_sum > 11 = 1, otherwise = 0

# SO1_2 - SO1_6: Intended behavior in each scenario. One analysis for each 
## scenario, where correct answers for all alternatives = 1, otherwise = 1.
## Scenario 1: if_else((sc1_1 == 1 & sc1_2 == 1 & sc1_3 == 1), 1, 0)
## One analysis for scenario 1 (SO1_2), one analysis for scenario 2 (SO1_3)...

# ------------------------------------
# SQ2. Do members of the public understand the main purpose of the advice 
## when experiencing symptoms of a newly onset respiratory tract infection 
## based on public health advice?

# SO2_1. Understanding of the advice. What is the aim of following the health
## advice?
## if_else(main_aim_advice == 'Protect people at higher risk of becoming severely ill', 1, 0)

# ------------------------------------
# SQ3. Does providing additional explanation of key terms used in infection 
### control advice facilitate understanding of those terms?

# SO3_1. Understanding relevant terms. Which of the following is included 
## under the term "infant"?
## if_else((infant_1 == TRUE & infant_2 == TRUE & infant_4 == TRUE), 1, 0)

# SO3_2. Understanding relevant terms. Which of the following is included 
## under the term "people at higher risk of becoming seriously ill"
## if_else((high_risk_1 == TRUE & high_risk_2 == TRUE & 
## high_risk_4 == TRUE & high_risk_5 == TRUE), 1, 0)

# -------------------------------------
# SQ4. Does explicitly adding a sentence stating that it is usually okay to go 
## to work or participate in activities if you feel up to it even with symptoms
## of a respiratory tract infection, improve understanding of, and correct 
## intended adherence to, infection control advice?
# SQ4 is interpreted based on analyses above

# -------------------------------------
# SQ5. Does explicitly adding a sentence stating that it is usually okay to go 
## to work or participate in activities if you feel up to it even with symptoms
## of a respiratory tract infection, and/or additional explanations of key 
## terms, affect trustworthiness and/or understandability of infection control
## advice?

# SO5_1. How trustworthy do you find the advice
## Comparing the arms based on an ordinale scale in variable: trustworthiness

# SO5_2. If I was looking for information about what to do when I had new 
## symptoms of respiratory illness, I would have found this advice useful.
## Comparing the arms based on an ordinale scale in variable: perceived_usefulness

# SO5_3. I would share this advice to a friend if they wanted to know what to 
## do in the case of new symptoms of respiratory illness.
## Comparing the arms based on an ordinale scale in variable: intention_to_share

# ------------------------------------------------------------------------------
# Start by running descriptive_analysis.R (to start a new result folder) and
## run prim_outcome.R (to ensure correct variables in panel)

# -------------------------------------------------------------------------
# Multiplicity across the secondary-outcome family
# -------------------------------------------------------------------------
## The protocol specifies 12 secondary analyses (across 6 research
## questions), each mirroring the primary analysis with 3 pairwise
## comparisons of active arms vs. V1_control -> 12 x 3 = 36 comparisons make
## up the secondary family. Define this once and reuse it consistently
## across all 12 analysis sections/scripts (see the multiplicity note after
## SO1_1 below for how to use it).
n_secondary_analyses       <- 12
n_comparisons_per_analysis <- 3
n_secondary_comparisons    <- n_secondary_analyses * n_comparisons_per_analysis  # 36

panel_sec <- panel |>
  mutate(sec_1 = if_else(scenario_sum > 11, 1, 0))

panel_sec |>
  group_by(arm) |>
  summarise(n = n(), n_meeting_threshold = sum(sec_1), pct = 100 * mean(sec_1))

# -------------------------------------------------------------------------
# SO1_1a. SATE -- log-binomial regression (RR of meeting the >=12/15
#         threshold, vs. V1_control)
# -------------------------------------------------------------------------
## Unlike the primary outcome (aggregated cbind(success, failure) counts per
## respondent), sec_1 is already a single 0/1 outcome per respondent, so this
## is a standard person-level log-binomial GLM. A Pearson-based dispersion
## check (as used for the primary outcome) is not informative here: a
## Bernoulli response cannot show over-dispersion in the way a binomial
## count can. HC3 robust SEs are used regardless, as protection against mild
## log-link misspecification (log-binomial models can be numerically fragile
## near the boundary of the parameter space).

fit_rr_so1_1 <- glm(
  sec_1 ~ arm,
  family = binomial(link = "log"),
  data = panel_sec
)

summary(fit_rr_so1_1)

V_so1_1 <- vcovHC(fit_rr_so1_1, type = "HC3")
coeftest(fit_rr_so1_1, vcov = V_so1_1)

# -------------------------------------------------------------------------
# SO1_1b. Bonferroni adjustment (within-analysis: 3 comparisons)
# -------------------------------------------------------------------------
## `m_so1_1` here is the number of comparisons THIS analysis contributes.
## To correct across the FULL secondary family instead, replace `m_so1_1`
## with `n_secondary_comparisons` when computing alpha_bonf_so1_1 /
## z_crit_so1_1 -- see the multiplicity note at the end of this section.

beta_so1_1 <- coef(fit_rr_so1_1)[-1]
se_so1_1   <- sqrt(diag(V_so1_1))[-1]

m_so1_1          <- length(beta_so1_1)   # 3 pairwise comparisons vs. control
alpha_bonf_so1_1 <- 0.05 / m_so1_1
z_crit_so1_1     <- qnorm(1 - alpha_bonf_so1_1 / 2)

z_value_so1_1 <- beta_so1_1 / se_so1_1
p_raw_so1_1   <- 2 * pnorm(abs(z_value_so1_1), lower.tail = FALSE)
p_bonf_so1_1  <- p.adjust(p_raw_so1_1, method = "bonferroni")

coefs_rr_so1_1 <- coef(fit_rr_so1_1)
rd_so1_1 <- purrr::map_dfr(seq_along(beta_so1_1) + 1, function(i) {
  out <- delta_rd(coefs_rr_so1_1, V_so1_1, i, ref_idx = 1)
  tibble(term = names(coefs_rr_so1_1)[i], RD = out["rd"], RD_SE = out["se"])
}) |>
  mutate(RD_CI_lower = RD - z_crit_so1_1 * RD_SE, RD_CI_upper = RD + z_crit_so1_1 * RD_SE)

resultat_so1_1 <- tibble(
  term      = names(beta_so1_1),
  RR        = exp(beta_so1_1),
  CI_lower  = exp(beta_so1_1 - z_crit_so1_1 * se_so1_1),
  CI_upper  = exp(beta_so1_1 + z_crit_so1_1 * se_so1_1),
  p_bonferroni = p_bonf_so1_1
) |>
  left_join(rd_so1_1, by = "term") |>
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

resultat_so1_1

# -------------------------------------------------------------------------
# SO1_1c. PATE -- population-weighted version (post-stratification weights
#         from weighting.R), mirroring prim_outcome.R section 6.
# -------------------------------------------------------------------------
## sec_1 already has one row per respondent, so (unlike the item-level PATE
## model for the primary outcome) there is no need to pivot to long format:
## svyglm() is fit directly, one row per participant.

panel_weighting <- readRDS("panel_weighting.rds") |>
  select(participant_id, weight)

panel_sec_w <- panel_sec |>
  left_join(panel_weighting, by = "participant_id")

design_w_so1_1 <- svydesign(ids = ~participant_id, weights = ~weight, data = panel_sec_w)

fit_pate_so1_1 <- svyglm(sec_1 ~ arm, design = design_w_so1_1, family = quasibinomial(link = "log"))

summary(fit_pate_so1_1)

beta_pate_so1_1 <- coef(fit_pate_so1_1)[-1]
se_pate_so1_1   <- sqrt(diag(vcov(fit_pate_so1_1)))[-1]

m_pate_so1_1          <- length(beta_pate_so1_1)
alpha_bonf_pate_so1_1 <- 0.05 / m_pate_so1_1
df_pate_so1_1         <- df.residual(fit_pate_so1_1)
t_crit_pate_so1_1     <- qt(1 - alpha_bonf_pate_so1_1 / 2, df = df_pate_so1_1)

t_value_pate_so1_1 <- beta_pate_so1_1 / se_pate_so1_1
p_raw_pate_so1_1   <- 2 * pt(abs(t_value_pate_so1_1), df = df_pate_so1_1, lower.tail = FALSE)
p_bonf_pate_so1_1  <- p.adjust(p_raw_pate_so1_1, method = "bonferroni")

coefs_pate_so1_1 <- coef(fit_pate_so1_1)
Vp_so1_1 <- vcov(fit_pate_so1_1)
rd_pate_so1_1 <- purrr::map_dfr(seq_along(beta_pate_so1_1) + 1, function(i) {
  out <- delta_rd(coefs_pate_so1_1, Vp_so1_1, i, ref_idx = 1)
  tibble(term = names(coefs_pate_so1_1)[i], RD = out["rd"], RD_SE = out["se"])
}) |>
  mutate(RD_CI_lower = RD - t_crit_pate_so1_1 * RD_SE, RD_CI_upper = RD + t_crit_pate_so1_1 * RD_SE)

resultat_pate_so1_1 <- tibble(
  term      = names(beta_pate_so1_1),
  RR        = exp(beta_pate_so1_1),
  CI_lower  = exp(beta_pate_so1_1 - t_crit_pate_so1_1 * se_pate_so1_1),
  CI_upper  = exp(beta_pate_so1_1 + t_crit_pate_so1_1 * se_pate_so1_1),
  p_bonferroni = p_bonf_pate_so1_1
) |>
  left_join(rd_pate_so1_1, by = "term") |>
  mutate(
    Comparison = recode(
      term,
      "armV2_sentence"             = "V2 (sentence) vs. V1 (control)",
      "armV3_definitions"          = "V3 (definitions) vs. V1 (control)",
      "armV4_sentence_definitions" = "V4 (sentence + definitions) vs. V1 (control)"
    ),
    `PATE RR (95% Bonferroni-adjusted CI)` = sprintf("%.2f (%.2f to %.2f)", RR, CI_lower, CI_upper),
    `PATE RD (95% Bonferroni-adjusted CI, pct. points)` = sprintf(
      "%.1f (%.1f to %.1f)", 100 * RD, 100 * RD_CI_lower, 100 * RD_CI_upper
    ),
    `Bonferroni-adjusted p value` = case_when(
      p_bonferroni < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_bonferroni)
    )
  ) |>
  select(Comparison, `PATE RR (95% Bonferroni-adjusted CI)`,
         `PATE RD (95% Bonferroni-adjusted CI, pct. points)`, `Bonferroni-adjusted p value`)

resultat_pate_so1_1

# -------------------------------------------------------------------------
# NOTE on multiplicity across all 12 secondary analyses
# -------------------------------------------------------------------------
## The Bonferroni correction above (m = 3) only accounts for the 3 pairwise
## comparisons WITHIN SO1_1, exactly as in prim_outcome.R -- it does NOT yet
## account for running 12 such analyses (36 comparisons in total). Three
## options, in increasing order of statistical power (decreasing
## conservativeness):
##
## 1. Strict Bonferroni across all 36 comparisons: replace `m_so1_1` /
##    `m_pate_so1_1` with `n_secondary_comparisons` (36) everywhere above.
##    Exact FWER control, but alpha = 0.05/36 = 0.0014 per test is very
##    conservative -- CIs widen and adjusted p-values inflate substantially.
##    Use if ALL 12 secondary analyses are pre-specified as confirmatory.
##
## 2. Holm-Bonferroni step-down across the pooled 36 raw p-values: also
##    controls the FWER exactly, but is uniformly more powerful than plain
##    Bonferroni. Requires pooling all 36 raw p-values FIRST, then adjusting
##    once (sketch below). A common compromise in practice: keep the CIs at
##    the Bonferroni (m = 36) width for a conservative interval estimate,
##    but report Holm-adjusted p-values for the significance decision.
##
## 3. Benjamini-Hochberg (FDR) across the pooled 36 raw p-values:
##    controls the expected proportion of false positives among results
##    declared "significant", rather than the probability of any false
##    positive at all. Noticeably more powerful than Bonferroni or Holm.
##    Appropriate if the 12 secondary analyses are treated as
##    hypothesis-generating/exploratory rather than confirmatory -- the
##    usual framing for secondary endpoints (ICH E9).
##
## Pooling sketch, to run once at the very end after all 12 sections'
## fit_rr_soX_Y / fit_pate_soX_Y models exist:
##
##   p_raw_all  <- c(p_raw_so1_1, p_raw_so1_2, ..., p_raw_so6_2)  # length 36
##   p_holm_all <- p.adjust(p_raw_all, method = "holm")
##   p_bh_all   <- p.adjust(p_raw_all, method = "BH")
##
## Decide (and pre-register in the SAP) which of the three is used BEFORE
## looking at results -- the choice materially changes which comparisons
## cross the significance threshold. Given 12 distinct research questions
## rather than 12 looks at the same question, Holm or BH across the pooled
## 36 tests is a defensible, less-conservative alternative to strict
## Bonferroni; confirm the choice with whoever owns the SAP.

# -------------------------------------------------------------------------
# SO1_2 - SO1_6. Intended behaviour within each individual scenario
# -------------------------------------------------------------------------
## One analysis per scenario (5 total): outcome = 1 if all 3 items for that
## scenario are answered correctly, else 0. Same modelling approach as
## SO1_1 (person-level log-binomial RR + Bonferroni, and the population-
## weighted PATE analogue), so this is refactored into one reusable
## function instead of repeating the ~130 lines from SO1_1 five times.

panel_sec <- panel_sec |>
  mutate(
    sec_2 = if_else(sc1_1 == 1 & sc1_2 == 1 & sc1_3 == 1, 1, 0),  # SO1_2: scenario 1
    sec_3 = if_else(sc2_1 == 1 & sc2_2 == 1 & sc2_3 == 1, 1, 0),  # SO1_3: scenario 2
    sec_4 = if_else(sc3_1 == 1 & sc3_2 == 1 & sc3_3 == 1, 1, 0),  # SO1_4: scenario 3
    sec_5 = if_else(sc4_1 == 1 & sc4_2 == 1 & sc4_3 == 1, 1, 0),  # SO1_5: scenario 4
    sec_6 = if_else(sc5_1 == 1 & sc5_2 == 1 & sc5_3 == 1, 1, 0)   # SO1_6: scenario 5
  )

panel_sec |>
  summarise(across(c(sec_2, sec_3, sec_4, sec_5, sec_6), ~ sum(.x, na.rm = TRUE)))

panel_sec |>
  group_by(arm) |>
  summarise(across(c(sec_2, sec_3, sec_4, sec_5, sec_6), mean), n = n())

## Reusable Comparison-label recoding, shared by every scenario analysis.
recode_comparison <- function(term) {
  recode(
    term,
    "armV2_sentence"             = "V2 (sentence) vs. V1 (control)",
    "armV3_definitions"          = "V3 (definitions) vs. V1 (control)",
    "armV4_sentence_definitions" = "V4 (sentence + definitions) vs. V1 (control)"
  )
}

## Fits both the SATE (log-binomial GLM, HC3 SEs) and PATE (svyglm,
## post-stratification-weighted) models for a single 0/1 outcome variable,
## and returns Bonferroni-adjusted RR/RD tables for both, mirroring SO1_1's
## logic exactly. `n_comparisons` sets the divisor for the Bonferroni
## correction: 3 for within-analysis correction (as used for SO1_1), or
## `n_secondary_comparisons` (36) to correct across the full secondary
## family instead -- see the multiplicity note above.
analyse_secondary_binary <- function(data, outcome, panel_weighting,
                                      n_comparisons = 3) {
  data <- data |> mutate(.outcome = .data[[outcome]])

  # ---- SATE: log-binomial GLM with HC3 robust SEs ----
  ## For high-prevalence ("ceiling") outcomes the log-binomial model can fail
  ## to converge (the log link requires all fitted probabilities <= 1, which
  ## becomes numerically fragile when the baseline risk is already close to
  ## 1). If that happens, fall back to a modified Poisson regression (Zou,
  ## 2004) with HC3 robust SEs: this has no upper-boundary constraint, gives
  ## the same RR interpretation, and remains valid (if slightly conservative)
  ## because the robust SE does not rely on the Poisson variance assumption.
  fit_rr <- tryCatch(
    glm(.outcome ~ arm, family = binomial(link = "log"), data = data),
    error = function(e) NULL
  )
  model_type <- "log-binomial"
  if (is.null(fit_rr) || !isTRUE(fit_rr$converged)) {
    message(sprintf(
      "Outcome '%s': log-binomial GLM did not converge (likely a high-prevalence/ceiling outcome) -- falling back to modified Poisson regression with HC3 robust SEs.",
      outcome
    ))
    fit_rr <- glm(.outcome ~ arm, family = poisson(link = "log"), data = data)
    model_type <- "Poisson (robust SE)"
  }
  V <- vcovHC(fit_rr, type = "HC3")

  beta <- coef(fit_rr)[-1]
  se   <- sqrt(diag(V))[-1]

  alpha_bonf <- 0.05 / n_comparisons
  z_crit     <- qnorm(1 - alpha_bonf / 2)

  z_value <- beta / se
  p_raw   <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
  p_bonf  <- p.adjust(p_raw, method = "bonferroni", n = n_comparisons)

  coefs <- coef(fit_rr)
  rd <- purrr::map_dfr(seq_along(beta) + 1, function(i) {
    out <- delta_rd(coefs, V, i, ref_idx = 1)
    tibble(term = names(coefs)[i], RD = out["rd"], RD_SE = out["se"])
  }) |>
    mutate(RD_CI_lower = RD - z_crit * RD_SE, RD_CI_upper = RD + z_crit * RD_SE)

  resultat <- tibble(
    term = names(beta), RR = exp(beta),
    CI_lower = exp(beta - z_crit * se), CI_upper = exp(beta + z_crit * se),
    p_raw = p_raw, p_bonferroni = p_bonf
  ) |>
    left_join(rd, by = "term") |>
    mutate(
      Comparison = recode_comparison(term),
      `RR (95% Bonferroni-adjusted CI)` = sprintf("%.2f (%.2f to %.2f)", RR, CI_lower, CI_upper),
      `RD (95% Bonferroni-adjusted CI, pct. points)` = sprintf(
        "%.1f (%.1f to %.1f)", 100 * RD, 100 * RD_CI_lower, 100 * RD_CI_upper
      ),
      `Bonferroni-adjusted p value` = case_when(
        p_bonferroni < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p_bonferroni)
      )
    ) |>
    select(term, Comparison, `RR (95% Bonferroni-adjusted CI)`,
           `RD (95% Bonferroni-adjusted CI, pct. points)`, `Bonferroni-adjusted p value`,
           p_raw)

  # ---- PATE: svyglm on post-stratification weights ----
  ## Same ceiling-outcome convergence risk as the SATE model above, so the
  ## same fallback applies: quasipoisson(link = "log") in place of
  ## quasibinomial(link = "log") when the latter fails to converge.
  data_w <- data |> left_join(panel_weighting, by = "participant_id")
  design_w <- svydesign(ids = ~participant_id, weights = ~weight, data = data_w)
  fit_pate <- tryCatch(
    svyglm(.outcome ~ arm, design = design_w, family = quasibinomial(link = "log")),
    error = function(e) NULL
  )
  model_type_pate <- "quasibinomial (log-binomial)"
  if (is.null(fit_pate) || !isTRUE(fit_pate$converged)) {
    message(sprintf(
      "Outcome '%s': PATE quasibinomial(log) svyglm did not converge -- falling back to quasipoisson(log).",
      outcome
    ))
    fit_pate <- svyglm(.outcome ~ arm, design = design_w, family = quasipoisson(link = "log"))
    model_type_pate <- "quasipoisson (robust SE)"
  }

  beta_pate <- coef(fit_pate)[-1]
  se_pate   <- sqrt(diag(vcov(fit_pate)))[-1]

  df_pate     <- df.residual(fit_pate)
  t_crit_pate <- qt(1 - alpha_bonf / 2, df = df_pate)

  t_value_pate <- beta_pate / se_pate
  p_raw_pate   <- 2 * pt(abs(t_value_pate), df = df_pate, lower.tail = FALSE)
  p_bonf_pate  <- p.adjust(p_raw_pate, method = "bonferroni", n = n_comparisons)

  coefs_pate <- coef(fit_pate)
  Vp <- vcov(fit_pate)
  rd_pate <- purrr::map_dfr(seq_along(beta_pate) + 1, function(i) {
    out <- delta_rd(coefs_pate, Vp, i, ref_idx = 1)
    tibble(term = names(coefs_pate)[i], RD = out["rd"], RD_SE = out["se"])
  }) |>
    mutate(RD_CI_lower = RD - t_crit_pate * RD_SE, RD_CI_upper = RD + t_crit_pate * RD_SE)

  resultat_pate <- tibble(
    term = names(beta_pate), RR = exp(beta_pate),
    CI_lower = exp(beta_pate - t_crit_pate * se_pate),
    CI_upper = exp(beta_pate + t_crit_pate * se_pate),
    p_raw = p_raw_pate, p_bonferroni = p_bonf_pate
  ) |>
    left_join(rd_pate, by = "term") |>
    mutate(
      Comparison = recode_comparison(term),
      `PATE RR (95% Bonferroni-adjusted CI)` = sprintf("%.2f (%.2f to %.2f)", RR, CI_lower, CI_upper),
      `PATE RD (95% Bonferroni-adjusted CI, pct. points)` = sprintf(
        "%.1f (%.1f to %.1f)", 100 * RD, 100 * RD_CI_lower, 100 * RD_CI_upper
      ),
      `Bonferroni-adjusted p value` = case_when(
        p_bonferroni < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p_bonferroni)
      )
    ) |>
    select(term, Comparison, `PATE RR (95% Bonferroni-adjusted CI)`,
           `PATE RD (95% Bonferroni-adjusted CI, pct. points)`, `Bonferroni-adjusted p value`,
           p_raw)

  list(fit_rr = fit_rr, fit_pate = fit_pate,
       model_type = model_type, model_type_pate = model_type_pate,
       resultat = resultat, resultat_pate = resultat_pate)
}

## Run SO1_2 - SO1_6 (scenarios 1-5) with a single Bonferroni divisor
## (n_comparisons) shared across the function calls. Kept at 3 (i.e.
## within-analysis correction, consistent with SO1_1 above) by default --
## change to n_secondary_comparisons to correct across all 12 analyses
## instead; see the multiplicity note above SO1_2.
scenario_outcomes <- c(
  SO1_2 = "sec_2",  # scenario 1
  SO1_3 = "sec_3",  # scenario 2
  SO1_4 = "sec_4",  # scenario 3
  SO1_5 = "sec_5",  # scenario 4
  SO1_6 = "sec_6"   # scenario 5
)

so1_2_6_results <- purrr::map(
  scenario_outcomes,
  ~ analyse_secondary_binary(panel_sec, .x, panel_weighting, n_comparisons = 3)
)

# SATE and PATE result tables per scenario, e.g.:
so1_2_6_results$SO1_2$resultat        # scenario 1, SATE
so1_2_6_results$SO1_2$resultat_pate   # scenario 1, PATE
so1_2_6_results$SO1_3$resultat        # scenario 2, SATE
so1_2_6_results$SO1_4$resultat        # scenario 3, SATE
so1_2_6_results$SO1_5$resultat        # scenario 4, SATE
so1_2_6_results$SO1_6$resultat        # scenario 5, SATE

# All 5 SATE tables stacked, with a Scenario column, if a single combined
# view/table is preferred over five separate ones:
so1_2_6_sate_combined <- purrr::imap_dfr(
  so1_2_6_results, ~ mutate(.x$resultat, Analysis = .y, .before = 1)
)
so1_2_6_sate_combined

so1_2_6_pate_combined <- purrr::imap_dfr(
  so1_2_6_results, ~ mutate(.x$resultat_pate, Analysis = .y, .before = 1)
)
so1_2_6_pate_combined

# -------------------------------------------------------------------------
# Pooled multiplicity adjustment across all 12 secondary analyses (sketch)
# -------------------------------------------------------------------------
## Once SO1_1 and SO1_2-SO1_6 (and the remaining SQ2-SQ6 analyses) have all
## been run, pool their raw p-values and apply Holm or BH ONCE across the
## full secondary family, per the multiplicity note above:
##
##   p_raw_all <- c(
##     p_raw_so1_1, p_raw_pate_so1_1,
##     purrr::map(so1_2_6_results, ~ .x$resultat$p_raw) |> unlist(),
##     purrr::map(so1_2_6_results, ~ .x$resultat_pate$p_raw) |> unlist()
##     # ... plus SQ2-SQ6 raw p-values once those scripts exist
##   )
##   p_holm_all <- p.adjust(p_raw_all, method = "holm")
##   p_bh_all   <- p.adjust(p_raw_all, method = "BH")
##
## (Now that all 12 secondary analyses are implemented, this pooling is
## carried out for real at the end of this script -- see the "Final
## combined secondary-outcome results" section.)

# -------------------------------------------------------------------------
# SQ2. Understanding of the main purpose of the advice
# -------------------------------------------------------------------------
## SO2_1. Outcome: correctly identifying the aim of the health advice
## (main_aim_advice == "Protect people at higher risk of becoming severely
## ill" = 1, otherwise 0). Single dichotomous item (not an all-3-correct
## composite as in SO1_2-SO1_6), but modelled the same way, so it reuses
## analyse_secondary_binary() defined above.

# Check for missingness in the source item before dichotomising -- any NA
# in main_aim_advice would silently become sec_7 = 0 via if_else() unless
# handled explicitly.
sum(is.na(panel_sec$main_aim_advice))

panel_sec <- panel_sec |>
  mutate(sec_7 = if_else(main_aim_advice == "Protect people at higher risk of becoming severely ill", 1, 0))

panel_sec |>
  group_by(arm) |>
  summarise(n = n(), n_correct = sum(sec_7), pct = 100 * mean(sec_7))

so2_1_results <- analyse_secondary_binary(panel_sec, "sec_7", panel_weighting, n_comparisons = 3)

resultat_so2_1      <- so2_1_results$resultat
resultat_pate_so2_1 <- so2_1_results$resultat_pate

resultat_so2_1
resultat_pate_so2_1

# -------------------------------------------------------------------------
# SQ3. Does providing additional explanation of key terms facilitate
#      understanding of those terms?
# -------------------------------------------------------------------------
## SO3_1 and SO3_2 are both "select all correct items" composites (like
## SO1_2-SO1_6), just applied to term-comprehension multi-select items
## (logical TRUE/FALSE columns) rather than scenario items. No missingness
## in the underlying infant_*/high_risk_* columns, so if_else() is safe to
## apply directly.

panel_sec <- panel_sec |>
  mutate(
    # SO3_1: term "infant" -- correct if infant_1, infant_2 and infant_4
    # (and only those) are identified as included
    sec_8 = if_else(infant_1 == TRUE & infant_2 == TRUE & infant_4 == TRUE, 1, 0),
    # SO3_2: term "people at higher risk of becoming seriously ill"
    sec_9 = if_else(
      high_risk_1 == TRUE & high_risk_2 == TRUE &
        high_risk_4 == TRUE & high_risk_5 == TRUE,
      1, 0
    )
  )

panel_sec |>
  group_by(arm) |>
  summarise(
    n = n(),
    n_correct_SO3_1 = sum(sec_8), pct_SO3_1 = 100 * mean(sec_8),
    n_correct_SO3_2 = sum(sec_9), pct_SO3_2 = 100 * mean(sec_9)
  )

so3_1_results <- analyse_secondary_binary(panel_sec, "sec_8", panel_weighting, n_comparisons = 3)
so3_2_results <- analyse_secondary_binary(panel_sec, "sec_9", panel_weighting, n_comparisons = 3)

resultat_so3_1      <- so3_1_results$resultat
resultat_pate_so3_1 <- so3_1_results$resultat_pate
resultat_so3_2      <- so3_2_results$resultat
resultat_pate_so3_2 <- so3_2_results$resultat_pate

resultat_so3_1
resultat_pate_so3_1
resultat_so3_2
resultat_pate_so3_2

# -------------------------------------------------------------------------
# SQ5. Trustworthiness and understandability of the advice
# -------------------------------------------------------------------------
## NB: the protocol text (see comments above, line 55-64) describes
## trustworthiness as a 1-10 Likert scale and perceived_usefulness /
## intention_to_share as 1-5 Likert scales. The actual variables in
## panel are ordered FACTORS with fewer levels: trustworthiness has 5
## levels ("Very little".."Very much"), and perceived_usefulness /
## intention_to_share each have 4 levels ("Strongly disagree".."Strongly
## agree") -- not numeric 1-10/1-5 scores. Flagging this discrepancy for the
## record.
##
## Decision (per discussion): dichotomise each into a "top-box" indicator
## and reuse the same log-binomial RR / Bonferroni framework as the other
## secondary analyses, for consistency and interpretability, at the cost of
## some information loss relative to modelling the full ordinal scale
## (e.g. a proportional-odds model would use all response categories).
## Top-box cut points:
##   SO5_1 trustworthiness:      "Much" or "Very much"          = 1, else 0
##   SO5_2 perceived_usefulness: "Agree" or "Strongly agree"     = 1, else 0
##   SO5_3 intention_to_share:   "Agree" or "Strongly agree"     = 1, else 0

sum(is.na(panel_sec$trustworthiness))
sum(is.na(panel_sec$perceived_usefulness))
sum(is.na(panel_sec$intention_to_share))

panel_sec <- panel_sec |>
  mutate(
    sec_10 = if_else(trustworthiness %in% c("Much", "Very much"), 1, 0),
    sec_11 = if_else(perceived_usefulness %in% c("Agree", "Strongly agree"), 1, 0),
    sec_12 = if_else(intention_to_share %in% c("Agree", "Strongly agree"), 1, 0)
  )

panel_sec |>
  group_by(arm) |>
  summarise(
    n = n(),
    pct_trust_topbox      = 100 * mean(sec_10),
    pct_useful_topbox     = 100 * mean(sec_11),
    pct_share_topbox      = 100 * mean(sec_12)
  )

so5_1_results <- analyse_secondary_binary(panel_sec, "sec_10", panel_weighting, n_comparisons = 3)
so5_2_results <- analyse_secondary_binary(panel_sec, "sec_11", panel_weighting, n_comparisons = 3)
so5_3_results <- analyse_secondary_binary(panel_sec, "sec_12", panel_weighting, n_comparisons = 3)

resultat_so5_1      <- so5_1_results$resultat
resultat_pate_so5_1 <- so5_1_results$resultat_pate
resultat_so5_2      <- so5_2_results$resultat
resultat_pate_so5_2 <- so5_2_results$resultat_pate
resultat_so5_3      <- so5_3_results$resultat
resultat_pate_so5_3 <- so5_3_results$resultat_pate

resultat_so5_1
resultat_pate_so5_1
resultat_so5_2
resultat_pate_so5_2
resultat_so5_3
resultat_pate_so5_3

# ===========================================================================
# FINAL COMBINED SECONDARY-OUTCOME RESULTS
# 12 analyses (36 pairwise comparisons); multiplicity-adjusted
# ===========================================================================
## All 12 secondary analyses are now implemented. This section:
##  1. Re-fits every analysis with the Bonferroni divisor set to
##     n_secondary_comparisons (36, i.e. the FULL secondary-outcome family,
##     not just the 3 comparisons within each individual analysis), giving
##     conservative Bonferroni-width 95% CIs for the PATE RR/RD.
##  2. Pools the 36 raw PATE p-values and applies Holm-Bonferroni step-down
##     correction ONCE across the pooled family -- more powerful than plain
##     Bonferroni while still controlling the family-wise error rate exactly
##     (see the multiplicity note earlier in this script for the rationale).
##  3. Assembles one combined table (mirroring Table 4 from prim_outcome.R:
##     crude/SATE n/N (%) per arm, calibrated/PATE % per arm, PATE RR/RD,
##     and an adjusted p-value) across all 36 comparisons, and exports it
##     as a publication-ready gt table.

n_secondary_analyses       <- 12
n_comparisons_per_analysis <- 3
n_secondary_comparisons    <- n_secondary_analyses * n_comparisons_per_analysis  # 36

secondary_outcomes_spec <- tibble::tribble(
  ~Analysis, ~outcome,  ~Description,
  "SO1_1",   "sec_1",   "SQ1: Meeting/exceeding threshold understanding (>=12/15 correct)",
  "SO1_2",   "sec_2",   "SQ1: Intended behaviour, scenario 1 (all 3 items correct)",
  "SO1_3",   "sec_3",   "SQ1: Intended behaviour, scenario 2 (all 3 items correct)",
  "SO1_4",   "sec_4",   "SQ1: Intended behaviour, scenario 3 (all 3 items correct)",
  "SO1_5",   "sec_5",   "SQ1: Intended behaviour, scenario 4 (all 3 items correct)",
  "SO1_6",   "sec_6",   "SQ1: Intended behaviour, scenario 5 (all 3 items correct)",
  "SO2_1",   "sec_7",   "SQ2: Correctly identifying the main aim of the advice",
  "SO3_1",   "sec_8",   "SQ3: Understanding the term \"infant\"",
  "SO3_2",   "sec_9",   "SQ3: Understanding the term \"people at higher risk\"",
  "SO5_1",   "sec_10",  "SQ5: Trustworthiness (top-box: Much/Very much)",
  "SO5_2",   "sec_11",  "SQ5: Perceived usefulness (top-box: Agree/Strongly agree)",
  "SO5_3",   "sec_12",  "SQ5: Intention to share (top-box: Agree/Strongly agree)"
)

stopifnot(nrow(secondary_outcomes_spec) == n_secondary_analyses)

# ---- 1. Re-fit all 12 analyses with the family-wide Bonferroni divisor ----
secondary_fits <- purrr::map(
  secondary_outcomes_spec$outcome,
  ~ analyse_secondary_binary(panel_sec, .x, panel_weighting, n_comparisons = n_secondary_comparisons)
) |>
  setNames(secondary_outcomes_spec$Analysis)

# Flag which analyses needed the modified-Poisson/quasipoisson fallback
# (ceiling outcomes), for the table footnote.
fallback_analyses <- secondary_outcomes_spec$Analysis[
  purrr::map_chr(secondary_fits, ~ .x$model_type_pate) != "quasibinomial (log-binomial)"
]
fallback_analyses

# ---- 2. Descriptive SATE n/N (%) and calibrated PATE % per arm ----
panel_sec_all_w <- panel_sec |> left_join(panel_weighting, by = "participant_id")
design_w_all <- svydesign(ids = ~participant_id, weights = ~weight, data = panel_sec_all_w)

describe_outcome <- function(outcome) {
  crude <- panel_sec |>
    group_by(arm) |>
    summarise(correct = sum(.data[[outcome]]), total = n(), .groups = "drop") |>
    mutate(sate_np = sprintf("%d/%d (%.1f%%)", correct, total, 100 * correct / total))

  calibrated <- svyby(as.formula(paste0("~", outcome)), ~arm, design_w_all, svymean) |>
    as_tibble() |>
    rename(pct = !!outcome) |>
    mutate(pate_pct = sprintf("%.1f%%", 100 * pct))

  list(crude = crude, calibrated = calibrated)
}

descriptives <- purrr::map(secondary_outcomes_spec$outcome, describe_outcome) |>
  setNames(secondary_outcomes_spec$Analysis)

term_to_arm <- c(
  "armV2_sentence"             = "V2_sentence",
  "armV3_definitions"          = "V3_definitions",
  "armV4_sentence_definitions" = "V4_sentence_definitions"
)

# ---- 3. Assemble the combined 36-row table ----
secondary_summary <- purrr::map2_dfr(
  secondary_outcomes_spec$Analysis, secondary_outcomes_spec$Description,
  function(analysis_label, description) {
    fit  <- secondary_fits[[analysis_label]]
    desc <- descriptives[[analysis_label]]
    control_sate_np  <- desc$crude$sate_np[desc$crude$arm == "V1_control"]
    control_pate_pct <- desc$calibrated$pate_pct[desc$calibrated$arm == "V1_control"]

    fit$resultat_pate |>
      mutate(
        comp_arm    = term_to_arm[term],
        Analysis    = analysis_label,
        Description = description,
        `SATE, n/N (%) control`          = control_sate_np,
        `SATE, n/N (%) intervention`     = desc$crude$sate_np[match(comp_arm, desc$crude$arm)],
        `PATE, calibrated % control`     = control_pate_pct,
        `PATE, calibrated % intervention` = desc$calibrated$pate_pct[match(comp_arm, desc$calibrated$arm)]
      )
  }
)

stopifnot(nrow(secondary_summary) == n_secondary_comparisons)

# Holm-Bonferroni: pool all 36 raw PATE p-values and adjust once. p.adjust()
# returns adjusted values aligned to the input order, so this can be added
# as a column directly.
secondary_summary <- secondary_summary |>
  mutate(
    p_holm = p.adjust(p_raw, method = "holm"),
    `Holm-Bonferroni-adjusted p value` = case_when(
      p_holm < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_holm)
    ),
    Group = paste0(Analysis, ". ", Description)
  ) |>
  select(
    Group, Comparison,
    `SATE, n/N (%) control`, `SATE, n/N (%) intervention`,
    `PATE, calibrated % control`, `PATE, calibrated % intervention`,
    `PATE RR (95% Bonferroni-width CI)` = `PATE RR (95% Bonferroni-adjusted CI)`,
    `PATE RD (95% Bonferroni-width CI, pct. points)` = `PATE RD (95% Bonferroni-adjusted CI, pct. points)`,
    `Holm-Bonferroni-adjusted p value`
  )

secondary_summary

# ---- 4. Publication-ready gt table (Table 5) ----
secondary_gt <- secondary_summary |>
  gt(groupname_col = "Group") |>
  tab_header(
    title = "Table 5. Secondary outcomes: population average treatment effect (PATE) vs. V1 (control)",
    subtitle = paste0(
      "12 secondary analyses (36 pairwise comparisons); post-stratification-weighted (age, gender, ",
      "education) log-binomial / modified-Poisson regression; 95% CIs at Bonferroni width ",
      "(alpha = 0.05/36); p values Holm-Bonferroni-adjusted across all 36 comparisons"
    )
  ) |>
  cols_label(
    Comparison = "Comparison",
    `SATE, n/N (%) control` = "SATE, n/N (%): control",
    `SATE, n/N (%) intervention` = "SATE, n/N (%): intervention",
    `PATE, calibrated % control` = "PATE, calibrated %: control",
    `PATE, calibrated % intervention` = "PATE, calibrated %: intervention",
    `PATE RR (95% Bonferroni-width CI)` = "PATE RR (95% Bonferroni-width CI)",
    `PATE RD (95% Bonferroni-width CI, pct. points)` = "PATE RD, pct. points (95% Bonferroni-width CI)",
    `Holm-Bonferroni-adjusted p value` = "Holm-Bonferroni-adjusted p value"
  ) |>
  cols_align(align = "left", columns = Comparison) |>
  cols_align(align = "center", columns = c(
    `SATE, n/N (%) control`, `SATE, n/N (%) intervention`,
    `PATE, calibrated % control`, `PATE, calibrated % intervention`,
    `PATE RR (95% Bonferroni-width CI)`,
    `PATE RD (95% Bonferroni-width CI, pct. points)`,
    `Holm-Bonferroni-adjusted p value`
  )) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) |>
  tab_footnote(
    footnote = paste(
      "SATE = sample average treatment effect (crude, unweighted trial counts).",
      "PATE = population average treatment effect: post-stratification weights",
      "from panel_weighting.rds (weighting.R), raked separately within each arm",
      "to national population margins for age group, gender, and education",
      "level (SSB). RD (risk difference) re-expresses the PATE RR on the",
      "absolute (percentage-point) scale via the delta method.",
      "CIs are reported at Bonferroni width (m = 36, the full secondary-outcome",
      "family) as a conservative interval estimate; p values are instead",
      "Holm-Bonferroni step-down adjusted across the same 36 comparisons,",
      "which is uniformly more powerful than Bonferroni while still",
      "controlling the family-wise error rate exactly.",
      if (length(fallback_analyses) > 0) {
        paste0(
          "For ", paste(fallback_analyses, collapse = ", "),
          ", the (quasi)binomial log-link model did not converge (ceiling/",
          "high-prevalence outcome), so a modified Poisson/quasipoisson",
          "regression with robust standard errors was used instead (Zou, 2004)."
        )
      } else {
        ""
      }
    )
  ) |>
  tab_options(
    table.font.size = px(10),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(11),
    column_labels.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.bottom.style = "solid"
  ) |>
  opt_table_font(font = "Times New Roman")

secondary_gt

gtsave_safe(secondary_gt, file.path(results_dir, "table5_secondary_outcomes.docx"))
gtsave_safe(secondary_gt, file.path(results_dir, "table5_secondary_outcomes.rtf"))
gtsave_safe(secondary_gt, file.path(results_dir, "table5_secondary_outcomes.html"))


