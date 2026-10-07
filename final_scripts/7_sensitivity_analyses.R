# Sensitivity analysis for the primary outcome: Manski (1990) worst-case /
# best-case bounds for missing outcome data.
## ------------------------------------------------------------------------
## The 140 dropouts (missing answer_time_ms -- no comprehension answers at
## all; see simulate_panel_test.R) are simply excluded from the primary
## analysis (prim_outcome.R). Manski bounds instead ask: across EVERY
## logically possible way those unobserved answers could have come out,
## how far could the estimate move? For a single arm's proportion correct,
## that means recoding all of that arm's missing items as wrong (the
## "worst case", a lower bound on the proportion) or as right (the "best
## case", an upper bound) -- the classic Manski (1990) bound.
##
## For a RATIO contrast between two arms (the RR vs. V1_control), the
## widest valid bound is NOT obtained by applying the same direction to
## both arms at once (that was what the previous version of this script
## did, fitting one model with everybody's missing data coded as failure
## and another with everybody's coded as success -- a valid-but-narrower
## sensitivity check, not the Manski bound for the RR itself). To bound a
## RATIO you must pair the extremes in OPPOSITE directions across the two
## arms being compared:
##   RR_lower = p_worst(intervention) / p_best(control)   (minimises RR)
##   RR_upper = p_best(intervention)  / p_worst(control)  (maximises RR)
## This is the traditional Manski-bound construction for a between-group
## contrast and is what's implemented below.
##
## Note these are deterministic bounds, not confidence intervals: they
## quantify uncertainty due to MISSINGNESS, not sampling error, so no SEs
## or p-values are attached.
## ------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(gt)

## Ensure results_dir exists (created by descreptive_analysis.R, line
## ~18-41) before running this script.

panel <- readRDS("panel_test.rds")

# Item-level correct/incorrect coding -- kept IDENTICAL to prim_outcome.R
# (including the two protocol deviations for scenario 3 & 4, item 3) so
# that the bounds below apply to the same primary outcome, just under
# different assumptions about the dropouts' unobserved answers.
panel <- panel |>
  mutate(
    sc1_1 = if_else(scenario1_item1 %in% c("Likely", "Very likely"), 1, 0),
    sc1_2 = if_else(scenario1_item2_careful %in% c("Likely", "Very likely"), 1, 0),
    sc1_3 = if_else(scenario1_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0),

    sc2_1 = if_else(scenario2_item1 %in% c("Likely", "Very likely"), 1, 0),
    sc2_2 = if_else(scenario2_item2_careful %in% c("Likely", "Very likely"), 1, 0),
    sc2_3 = if_else(scenario2_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0),

    sc3_1 = if_else(scenario3_item1 %in% c("Very unlikely", "Unlikely"), 1, 0),
    sc3_2 = if_else(scenario3_item2_careful %in% c("Likely", "Very likely"), 1, 0),
    # OBS! Deviation from the protocol: 4 & 5 (not 1 & 2) -- see prim_outcome.R
    sc3_3 = if_else(scenario3_item3_normal %in% c("Likely", "Very likely"), 1, 0),

    sc4_1 = if_else(scenario4_item1 %in% c("Very unlikely", "Unlikely"), 1, 0),
    sc4_2 = if_else(scenario4_item2_careful %in% c("Likely", "Very likely"), 1, 0),
    # OBS! Deviation from the protocol: 4 & 5 (not 1 & 2) -- see prim_outcome.R
    sc4_3 = if_else(scenario4_item3_normal %in% c("Likely", "Very likely"), 1, 0),

    sc5_1 = if_else(scenario5_item1 %in% c("Very unlikely", "Unlikely"), 1, 0),
    sc5_2 = if_else(scenario5_item2_careful %in% c("Likely", "Very likely"), 1, 0),
    sc5_3 = if_else(scenario5_item3_normal %in% c("Very unlikely", "Unlikely"), 1, 0),

    arm             = relevel(factor(arm), ref = "V1_control"),
    missing_outcome = is.na(answer_time_ms)
  )

sc_cols <- c(
  "sc1_1", "sc1_2", "sc1_3", "sc2_1", "sc2_2", "sc2_3",
  "sc3_1", "sc3_2", "sc3_3", "sc4_1", "sc4_2", "sc4_3",
  "sc5_1", "sc5_2", "sc5_3"
)
n_items <- length(sc_cols)   # 15

# -------------------------------------------------------------------------
# 1. Per-respondent worst-case / best-case item counts
# -------------------------------------------------------------------------
# For dropouts every sc_* column is NA, so the observed sum is computed
# with na.rm = TRUE (giving 0 for a fully-missing row) and then explicitly
# overridden to the two extremes for that row, rather than letting NAs
# silently fall through.
panel <- panel |>
  mutate(
    scenario_sum_observed = rowSums(across(all_of(sc_cols)), na.rm = TRUE),
    scenario_sum_worst     = if_else(missing_outcome, 0L, scenario_sum_observed),
    scenario_sum_best      = if_else(missing_outcome, as.integer(n_items), scenario_sum_observed)
  )

# -------------------------------------------------------------------------
# 2. Per-arm Manski bounds on the proportion of correct answers
#    p_worst = assume every dropout's 15 answers are ALL wrong (lower bound)
#    p_best  = assume every dropout's 15 answers are ALL right (upper bound)
# -------------------------------------------------------------------------
arm_bounds <- panel |>
  group_by(arm) |>
  summarise(
    N         = n(),
    n_missing = sum(missing_outcome),
    sum_worst = sum(scenario_sum_worst),
    sum_best  = sum(scenario_sum_best),
    .groups   = "drop"
  ) |>
  mutate(
    pct_missing = 100 * n_missing / N,
    p_worst     = sum_worst / (N * n_items),
    p_best      = sum_best  / (N * n_items)
  )

arm_bounds

# -------------------------------------------------------------------------
# 3. Traditional Manski bounds for the RR and RD vs. V1_control
#    (opposite-direction pairing across arms -- see header note)
# -------------------------------------------------------------------------
control <- arm_bounds |> filter(arm == "V1_control")

manski <- arm_bounds |>
  filter(arm != "V1_control") |>
  mutate(
    RR_lower = p_worst / control$p_best,
    RR_upper = p_best  / control$p_worst,
    RD_lower = p_worst - control$p_best,
    RD_upper = p_best  - control$p_worst,
    Comparison = recode(
      as.character(arm),
      "V2_sentence"             = "V2 (sentence) vs. V1 (control)",
      "V3_definitions"          = "V3 (definitions) vs. V1 (control)",
      "V4_sentence_definitions" = "V4 (sentence + definitions) vs. V1 (control)"
    )
  )

resultat_manski <- manski |>
  mutate(
    `RR Manski bounds`              = sprintf("%.2f to %.2f", RR_lower, RR_upper),
    `RD Manski bounds, pct. points` = sprintf("%.1f to %.1f", 100 * RD_lower, 100 * RD_upper),
    `Includes RR = 1?`              = if_else(RR_lower <= 1 & RR_upper >= 1, "Yes", "No")
  ) |>
  select(Comparison, `RR Manski bounds`, `RD Manski bounds, pct. points`, `Includes RR = 1?`)

resultat_manski

# -------------------------------------------------------------------------
# 4. Publication-ready table (gt)
# -------------------------------------------------------------------------
gtsave_safe <- function(data, path) {
  tmp <- tempfile(fileext = paste0(".", tools::file_ext(path)))
  gtsave(data, tmp)
  file.copy(tmp, path, overwrite = TRUE)
  file.remove(tmp)
  invisible(path)
}

manski_gt <- resultat_manski |>
  gt() |>
  tab_header(
    title    = "Table 7. Manski worst-case/best-case bounds for the primary outcome",
    subtitle = sprintf(
      "Bounds on the RR/RD of a correct comprehension answer vs. V1 (control), allowing the %d dropouts' unobserved answers to take any value",
      sum(arm_bounds$n_missing)
    )
  ) |>
  cols_label(
    Comparison = "Comparison",
    `RR Manski bounds` = "RR Manski bounds (worst case to best case)",
    `RD Manski bounds, pct. points` = "RD Manski bounds, pct. points (worst case to best case)",
    `Includes RR = 1?` = "Bounds include RR = 1?"
  ) |>
  cols_align(align = "left", columns = Comparison) |>
  cols_align(align = "center", columns = c(
    `RR Manski bounds`, `RD Manski bounds, pct. points`, `Includes RR = 1?`
  )) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_footnote(
    footnote = paste(
      "Manski (1990) bounds: deterministic worst-case/best-case limits on the",
      "estimate under complete uncertainty about missing outcomes, not",
      "statistical confidence intervals -- no SE or p-value is attached.",
      "Worst case = all of an arm's missing items scored wrong; best case =",
      "all scored right. RR/RD bounds pair the two arms' extremes in",
      "opposite directions (the arm's own worst case against the other",
      "arm's best case, and vice versa) to obtain the widest valid bound",
      "on the between-arm contrast."
    )
  ) |>
  tab_footnote(
    footnote = sprintf(
      paste(
        "Missing data (dropouts, no comprehension answers at all): %.1f%%",
        "(V1), %.1f%% (V2), %.1f%% (V3), %.1f%% (V4)."
      ),
      arm_bounds$pct_missing[arm_bounds$arm == "V1_control"],
      arm_bounds$pct_missing[arm_bounds$arm == "V2_sentence"],
      arm_bounds$pct_missing[arm_bounds$arm == "V3_definitions"],
      arm_bounds$pct_missing[arm_bounds$arm == "V4_sentence_definitions"]
    )
  ) |>
  tab_options(
    table.font.size            = px(12),
    heading.title.font.size    = px(14),
    heading.subtitle.font.size = px(12),
    column_labels.font.weight  = "bold",
    table.border.top.style     = "solid",
    table.border.bottom.style  = "solid"
  ) |>
  opt_table_font(font = "Times New Roman")

manski_gt

gtsave_safe(manski_gt, file.path(results_dir, "table7_manski_bounds.docx"))
gtsave_safe(manski_gt, file.path(results_dir, "table7_manski_bounds.rtf"))
gtsave_safe(manski_gt, file.path(results_dir, "table7_manski_bounds.html"))

# -------------------------------------------------------------------------
# 5. Figure: Manski RR bounds vs. RR = 1
# -------------------------------------------------------------------------
manski_plot_df <- manski |>
  mutate(Comparison = factor(Comparison, levels = rev(Comparison)))

manski_plot <- ggplot(manski_plot_df, aes(y = Comparison)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_segment(aes(x = RR_lower, xend = RR_upper, yend = Comparison), linewidth = 1.1) +
  geom_point(aes(x = RR_lower), shape = 17, size = 2.8) +
  geom_point(aes(x = RR_upper), shape = 17, size = 2.8) +
  geom_text(
    aes(x = RR_lower, label = sprintf("%.2f", RR_lower)),
    vjust = -1, size = 3.2
  ) +
  geom_text(
    aes(x = RR_upper, label = sprintf("%.2f", RR_upper)),
    vjust = -1, size = 3.2
  ) +
  scale_x_log10() +
  labs(
    x = "Relative risk, Manski worst case to best case (log scale)",
    y = NULL,
    title = "Manski bounds for the primary-outcome RR vs. V1 (control)",
    subtitle = sprintf(
      "Allowing all %d dropouts' unobserved comprehension answers to take any value",
      sum(arm_bounds$n_missing)
    )
  ) +
  theme_minimal()

manski_plot

ggsave(
  file.path(results_dir, "figure_manski_bounds.png"),
  manski_plot, width = 8, height = 4.5, dpi = 150, bg = "white"
)

#-------------------------------------------------------------------------------
