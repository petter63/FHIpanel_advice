# Primary outcome, all reference arms
## Same analysis as 4_prim_outcome.R, section 1-5 (SATE: relative risk of a
## correct comprehension answer out of 15, via log-binomial regression),
## repeated once per arm as the reference/control arm the other three are
## compared against. Produces 4 versions of Table 5, one per reference arm
## -- so every pairwise comparison between arms is available, not just the
## 3 comparisons vs. V1_control that 4_prim_outcome.R's Table 3 covers.
## The PATE (population-weighted) analysis (4_prim_outcome.R section 6
## onwards) is NOT repeated here -- this script only reruns the primary
## SATE RR/RD model.
##
## NOTE on multiplicity: each table below is still only Bonferroni-adjusted
## for the 3 comparisons vs. THAT table's reference arm, exactly as in
## 4_prim_outcome.R. Looked at together, the 4 tables contain 12 rows but
## only 6 unique arm pairs (each pair is estimated twice, once from each
## arm's side as reference) -- if we interpret results ACROSS the 4
## tables rather than reading one at a time, the effective number of
## comparisons under consideration is larger than 3 and a stricter
## correction (e.g. Bonferroni for the 6 unique pairs, ignoring the
## mirrored duplicates) should be considered instead.

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(sandwich)
library(lmtest)
library(gt)

# -------------------------------------------------------------------------
## Ensure to create a new folder and log by running 2_descreptive_analysis.R
## before this script (same prerequisite as 4_prim_outcome.R, for results_dir).

# -------------------------------------------------------------------------
# 1. Code correct answers and build the outcome (successes / failures out of 15)
#    -- identical to 4_prim_outcome.R section 1, but `arm` is left as plain
#    factor here (not releveled yet); each reference arm is set inside the
#    loop below instead.
# -------------------------------------------------------------------------

panel_base <- readRDS("panel_data.rds") |>
  filter(!is.na(answer_time_ms)) |>
  mutate(sc1_1 = if_else(scenario1_item1 %in% c("Likely", "Very likely"), 1, 0),
         sc1_2 = if_else(scenario1_item2_careful %in% c("Likely", "Very likely"), 1, 0),
         sc1_3 = if_else(scenario1_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0),

         sc2_1 = if_else(scenario2_item1 %in% c("Likely", "Very likely"), 1, 0),
         sc2_2 = if_else(scenario2_item2_careful %in% c("Likely", "Very likely"), 1, 0),
         sc2_3 = if_else(scenario2_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0),

         sc3_1 = if_else(scenario3_item1 %in% c("Very unlikely", "Unlikely"), 1, 0),
         sc3_2 = if_else(scenario3_item2_careful %in% c("Likely", "Very likely"), 1, 0),

         # OBS! Deviation from the protocol: 4 & 5 (not 1 & 2)
         sc3_3 = if_else(scenario3_item3_normal %in% c("Likely", "Very likely"), 1, 0),

         sc4_1 = if_else(scenario4_item1 %in% c("Very unlikely", "Unlikely"), 1, 0),
         sc4_2 = if_else(scenario4_item2_careful %in% c("Likely", "Very likely"), 1, 0),

         # OBS! Deviation from the protocol: 4 & 5 (not 1 & 2)
         sc4_3 = if_else(scenario4_item3_normal %in% c("Likely", "Very likely"), 1, 0),

         sc5_1 = if_else(scenario5_item1 %in% c("Very unlikely", "Unlikely"), 1, 0),
         sc5_2 = if_else(scenario5_item2_careful %in% c("Likely", "Very likely"), 1, 0),
         sc5_3 = if_else(scenario5_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0)) |>

  # Sum of all 15 items
  mutate(scenario_sum = rowSums(across(c(sc1_1, sc1_2, sc1_3,
                                          sc2_1, sc2_2, sc2_3,
                                          sc3_1, sc3_2, sc3_3,
                                          sc4_1, sc4_2, sc4_3,
                                          sc5_1, sc5_2, sc5_3)))) |>
  mutate(scenario_fail = 15 - scenario_sum)

range(panel_base$scenario_sum)
range(panel_base$scenario_fail)

arm_levels <- c("V1_control", "V2_sentence", "V3_definitions", "V4_sentence_definitions")
arm_full_label <- c(
  V1_control              = "V1 (control)",
  V2_sentence              = "V2 (sentence)",
  V3_definitions           = "V3 (definitions)",
  V4_sentence_definitions  = "V4 (sentence + definitions)"
)

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

## Risk difference (RD), re-expressed from the RR -- identical to
## 4_prim_outcome.R's delta_rd(): RD is a non-linear function of two
## correlated model parameters (control-arm log-risk b0, log-RR b1), so a
## delta method propagates uncertainty in BOTH into the SE of RD.
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

# -------------------------------------------------------------------------
# 2-5 wrapped into one function: log-binomial RR, over-dispersion check,
# Bonferroni-adjusted CIs/RD, and the publication-ready Table 5 -- run once
# per reference arm.
# -------------------------------------------------------------------------
run_primary_vs_ref <- function(panel_base, ref_arm) {
  panel <- panel_base |> mutate(arm = relevel(factor(arm), ref = ref_arm))

  fit_rr <- glm(
    cbind(scenario_sum, scenario_fail) ~ arm,
    family = binomial(link = "log"),
    data = panel
  )

  dispersion <- sum(residuals(fit_rr, type = "pearson")^2) / df.residual(fit_rr)
  use_robust_se <- dispersion > 1.2
  V <- if (use_robust_se) vcovHC(fit_rr, type = "HC3") else vcov(fit_rr)

  beta <- coef(fit_rr)[-1]
  se   <- sqrt(diag(V))[-1]

  m <- length(beta)                 # number of comparisons vs. this reference (3)
  alpha_bonf <- 0.05 / m
  z_crit <- qnorm(1 - alpha_bonf / 2)

  z_value <- beta / se
  p_raw   <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
  p_bonf  <- p.adjust(p_raw, method = "bonferroni")

  coefs_rr <- coef(fit_rr)
  rd_sate <- purrr::map_dfr(seq_along(beta) + 1, function(i) {
    out <- delta_rd(coefs_rr, V, i, ref_idx = 1)
    tibble(term = names(coefs_rr)[i], RD = out["rd"], RD_SE = out["se"])
  }) |>
    mutate(RD_CI_lower = RD - z_crit * RD_SE, RD_CI_upper = RD + z_crit * RD_SE)

  comp_arms   <- setdiff(arm_levels, ref_arm)
  term_labels <- setNames(
    paste(arm_full_label[comp_arms], "vs.", arm_full_label[ref_arm]),
    paste0("arm", comp_arms)
  )

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
      Comparison = term_labels[term],
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

  table5_gt <- resultat |>
    gt() |>
    tab_header(
      title = "Table 5. Relative risk of a correct comprehension answer",
      subtitle = sprintf(
        "Log-binomial regression vs. %s (reference); %s; Bonferroni-adjusted for %d comparisons",
        arm_full_label[ref_arm],
        if (use_robust_se) "robust (HC3) SEs used due to over-dispersion" else "model-based SEs",
        m
      )
    ) |>
    cols_label(
      Comparison = "Comparison",
      `RR (95% Bonferroni-adjusted CI)` = "RR (95% Bonferroni-adjusted CI)",
      `RD (95% Bonferroni-adjusted CI, pct. points)` = "RD, pct. points (95% Bonferroni-adjusted CI)",
      `Bonferroni-adjusted p value` = "Bonferroni-adjusted p value"
    ) |>
    cols_align(align = "left", columns = Comparison) |>
    cols_align(align = "center", columns = c(`RR (95% Bonferroni-adjusted CI)`,
                                              `RD (95% Bonferroni-adjusted CI, pct. points)`,
                                              `Bonferroni-adjusted p value`)) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) |>
    tab_footnote(
      footnote = sprintf(
        paste(
          "Outcome: number of correct answers out of 15 comprehension items.",
          "Reference (control) arm for this table: %s. Dispersion statistic = %.2f.",
          "RD (risk difference) is the RR re-expressed on the absolute",
          "(percentage-point) scale via the delta method, propagating",
          "uncertainty in both the reference-arm risk and the RR."
        ),
        arm_full_label[ref_arm], dispersion
      )
    ) |>
    tab_footnote(
      footnote = paste(
        "This table is one of 4 reference-arm versions of Table 5 (see",
        "4b_prim_outcome_allref.R); each is Bonferroni-adjusted only for its",
        "own 3 comparisons -- reading multiple versions together implies",
        "more comparisons than any single table's correction accounts for."
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

  list(
    ref_arm       = ref_arm,
    fit_rr        = fit_rr,
    dispersion    = dispersion,
    use_robust_se = use_robust_se,
    resultat      = resultat,
    table5_gt     = table5_gt
  )
}

# -------------------------------------------------------------------------
# Run all 4 reference-arm versions and save each as its own Table 5
# -------------------------------------------------------------------------
results_allref <- map(arm_levels, ~ run_primary_vs_ref(panel_base, .x))
names(results_allref) <- arm_levels

iwalk(results_allref, function(res, ref_arm) {
  print(res$table5_gt)

  file_stub <- sprintf("table5_primary_outcome_RR_ref_%s", ref_arm)
  gtsave_safe(res$table5_gt, file.path(results_dir, paste0(file_stub, ".docx")))
  gtsave_safe(res$table5_gt, file.path(results_dir, paste0(file_stub, ".rtf")))
  gtsave_safe(res$table5_gt, file.path(results_dir, paste0(file_stub, ".html")))
})

# Quick overview of all 12 rows (4 reference arms x 3 comparisons each)
# across the 4 tables, for a single at-a-glance sanity check.
map_dfr(results_allref, "resultat", .id = "reference_arm")

#-------------------------------------------------------------------------------
