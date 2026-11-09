#-------------------------------------------------------------------------------
# Meta-analysis of study 1 and study 2: primary outcome
#-------------------------------------------------------------------------------
# Pools the effect of V2, V3 and V4 vs. V1 (control) on a correct comprehension
# answer across the two studies, using the estimates saved by
# study1_4_prim_outcome.R and study2_4_prim_outcome.R:
#   - SA: sample-average log RR (log-binomial regression, unweighted)
#   - PA: population-average log RR (post-stratification weighted)
#   - RD (risk difference) via the delta method, pooled on the absolute scale
#
# Method: inverse-variance weighted meta-analysis of log RR (and RD).
#   - Fixed-effect pooling only. A random-effects model is not used: it needs
#     at least three studies (two parameters, mean and between-study variance,
#     are estimated), and only two studies are available.
#   - Heterogeneity: Cochran's Q, I^2 (descriptive only with k = 2).
#   - Output: the table (Table 8) is the reported result. The forest plot is
#     for internal use only.
#   - Pooled CIs and p values use the same Bonferroni adjustment as the
#     single-study analyses (m = 3 comparisons vs. V1).
#
# Assumptions to be aware of:
#   - V1-V4 are the same interventions in both studies, and arm labels refer to
#     the real (unblinded) allocation in both. 
#   - The two studies are independent (separate participants).
#   - Study 1 and study 2 differ in inclusion handling (see the organize-data
#     scripts); the estimates pooled here are from analysed participants only.
#
# Run study1_1 to study1_4 and study2_1 to study2_4 first.
#-------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(gt)

# Set to TRUE once the study 2 allocation has been unblinded and the real
# file -> arm mapping is used in study2_1_organize_data.R.
study2_unblinded <- FALSE

if (!study2_unblinded) {
  warning(
    "Study 2 allocation is still BLINDED (random arms). The pooled estimates ",
    "below are NOT interpretable until study 2 is unblinded.",
    call. = FALSE
  )
}

run_id      <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- file.path("results", "meta", run_id)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

gtsave_safe <- function(data, path) {
  tmp <- tempfile(fileext = paste0(".", tools::file_ext(path)))
  gtsave(data, tmp)
  file.copy(tmp, path, overwrite = TRUE)
  file.remove(tmp)
  invisible(path)
}

# -------------------------------------------------------------------------
# 1. Read per-study estimates
# -------------------------------------------------------------------------
est <- bind_rows(
  mutate(readRDS("study1_primary_estimates.rds"), study = "Study 1"),
  mutate(readRDS("study2_primary_estimates.rds"), study = "Study 2")
)

stopifnot(
  setequal(unique(est$term[est$study == "Study 1"]), unique(est$term[est$study == "Study 2"])),
  all(is.finite(est$logRR)), all(is.finite(est$SE)), all(est$SE > 0)
)

m_comparisons <- n_distinct(est$term)
z_crit <- qnorm(1 - (0.05 / m_comparisons) / 2)

comparison_labels <- c(
  armV2_sentence             = "V2 (sentence) vs. V1 (control)",
  armV3_definitions          = "V3 (definitions) vs. V1 (control)",
  armV4_sentence_definitions = "V4 (sentence + definitions) vs. V1 (control)"
)

# -------------------------------------------------------------------------
# 2. Inverse-variance pooling
# -------------------------------------------------------------------------
pool_ivw <- function(y, se, m) {
  w   <- 1 / se^2
  k   <- length(y)
  fe  <- sum(w * y) / sum(w)
  se_fe <- sqrt(1 / sum(w))

  Q   <- sum(w * (y - fe)^2)
  df  <- k - 1
  I2  <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0

  p_fe <- 2 * pnorm(abs(fe / se_fe), lower.tail = FALSE)

  tibble(
    est_fe = fe, se_fe = se_fe, p_fe_bonf = min(1, p_fe * m),
    Q = Q, p_Q = pchisq(Q, df, lower.tail = FALSE), I2 = I2
  )
}

pooled_rr <- est |>
  group_by(estimand, term) |>
  group_modify(~ pool_ivw(.x$logRR, .x$SE, m_comparisons)) |>
  ungroup()

pooled_rd <- est |>
  group_by(estimand, term) |>
  group_modify(~ pool_ivw(.x$RD, .x$RD_SE, m_comparisons)) |>
  ungroup()

# -------------------------------------------------------------------------
# 3. Results table
# -------------------------------------------------------------------------
fmt_rr <- function(logrr, se) {
  sprintf("%.2f (%.2f to %.2f)", exp(logrr), exp(logrr - z_crit * se), exp(logrr + z_crit * se))
}
fmt_rd <- function(rd, se) {
  sprintf("%.1f (%.1f to %.1f)", 100 * rd, 100 * (rd - z_crit * se), 100 * (rd + z_crit * se))
}
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

study_cols <- est |>
  mutate(rr = fmt_rr(logRR, SE)) |>
  select(estimand, term, study, rr) |>
  pivot_wider(names_from = study, values_from = rr)

meta_table <- pooled_rr |>
  left_join(
    pooled_rd |> transmute(estimand, term, rd_fe = est_fe, rd_se_fe = se_fe),
    by = c("estimand", "term")
  ) |>
  left_join(study_cols, by = c("estimand", "term")) |>
  mutate(
    Estimand   = if_else(estimand == "SA", "Sample average (SA)", "Population average (PA)"),
    Comparison = comparison_labels[term],
    `Pooled RR, fixed effect (95% Bonferroni-adjusted CI)` = fmt_rr(est_fe, se_fe),
    `Pooled RD, pct. points (95% Bonferroni-adjusted CI)`  = fmt_rd(rd_fe, rd_se_fe),
    `Bonferroni-adjusted p value`                           = fmt_p(p_fe_bonf),
    `Q (p)`  = sprintf("%.2f (%.2f)", Q, p_Q),
    `I2 (%)` = sprintf("%.0f", I2)
  ) |>
  arrange(desc(estimand), term) |>
  select(
    Estimand, Comparison,
    `Study 1 RR (95% CI)` = `Study 1`, `Study 2 RR (95% CI)` = `Study 2`,
    `Pooled RR, fixed effect (95% Bonferroni-adjusted CI)`,
    `Pooled RD, pct. points (95% Bonferroni-adjusted CI)`,
    `Bonferroni-adjusted p value`, `Q (p)`, `I2 (%)`
  )

meta_table

meta_gt <- meta_table |>
  gt(groupname_col = "Estimand") |>
  tab_header(
    title = "Table 8. Meta-analysis of study 1 and study 2: effect on a correct comprehension answer",
    subtitle = sprintf(
      "Inverse-variance fixed-effect pooling of log RR vs. V1 (control); Bonferroni-adjusted for %d comparisons",
      m_comparisons
    )
  ) |>
  cols_align(align = "left", columns = Comparison) |>
  cols_align(align = "center", columns = -Comparison) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) |>
  tab_footnote(
    footnote = paste(
      "SA = sample average (unweighted); PA = population average (post-stratification weighted).",
      "Study-specific RRs are the unadjusted estimates from study1_4_prim_outcome.R and",
      "study2_4_prim_outcome.R. RD is pooled on the absolute scale using each study's",
      "delta-method SE.",
      "Q = Cochran's Q, I2 = share of variability due to between-study heterogeneity;",
      "both are uninformative with only two studies.",
      "Fixed-effect pooling only; a random-effects model requires at least three studies."
    )
  ) |>
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    column_labels.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.bottom.style = "solid"
  ) |>
  opt_table_font(font = "Times New Roman")

meta_gt

gtsave_safe(meta_gt, file.path(results_dir, "table8_meta_analysis.docx"))
gtsave_safe(meta_gt, file.path(results_dir, "table8_meta_analysis.rtf"))
gtsave_safe(meta_gt, file.path(results_dir, "table8_meta_analysis.html"))

# -------------------------------------------------------------------------
# 4. Forest plot (internal use only -- not saved; Table 8 is the reported result)
# -------------------------------------------------------------------------
row_levels <- c("Pooled (fixed effect)", "Study 2", "Study 1")

forest_df <- bind_rows(
  est |>
    transmute(estimand, term, row = study, type = "Study",
              rr = exp(logRR), lower = exp(logRR - z_crit * SE), upper = exp(logRR + z_crit * SE)),
  pooled_rr |>
    transmute(estimand, term, row = "Pooled (fixed effect)", type = "Pooled",
              rr = exp(est_fe), lower = exp(est_fe - z_crit * se_fe), upper = exp(est_fe + z_crit * se_fe))
) |>
  mutate(
    row        = factor(row, levels = row_levels),
    Comparison = factor(comparison_labels[term], levels = comparison_labels),
    estimand   = factor(estimand, levels = c("SA", "PA"),
                        labels = c("Sample average (SA)", "Population average (PA)"))
  )

forest_plot <- ggplot(forest_df, aes(x = rr, y = row, shape = type)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.5) +
  scale_shape_manual(values = c(Study = 16, Pooled = 18), guide = "none") +
  scale_x_log10() +
  facet_grid(Comparison ~ estimand, labeller = label_wrap_gen(width = 24)) +
  labs(
    x = "Relative risk of a correct answer (log scale), 95% Bonferroni-adjusted CI",
    y = NULL,
    title = "Meta-analysis of study 1 and study 2",
    subtitle = "Primary outcome, vs. V1 (control)"
  ) +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0))

forest_plot
