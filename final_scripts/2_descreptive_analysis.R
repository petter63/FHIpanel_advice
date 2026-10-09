#-------------------------------------------------------------------------------
# Project: FHI-panel effectiveness of infection prevention recommendations
#-------------------------------------------------------------------------------
# Load and describe study data
# NB: the data frame is named "panel" from here on. 
library(skimr)

panel <- readRDS("panel_data.rds")
str(panel)
skimr::skim(panel)

#-------------------------------------------------------------------------------
# Generating result subfolder for each run
#-------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)

# Each run gets its own timestamped subfolder under results/, so output from
# different runs is kept separate and document when a given result was produced.
run_id      <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- file.path("results", run_id)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Record run metadata (timestamp, R version, git commit) for traceability.
git_commit <- tryCatch(
  system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE),
  error = function(e) NA_character_
)
if (length(git_commit) == 0) git_commit <- NA_character_

writeLines(
  c(
    paste("Run timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("R version:", R.version.string),
    paste("Git commit:", git_commit),
    paste("Script:", "descreptive_analysis.R")
  ),
  file.path(results_dir, "run_info.txt")
)

n_total <- nrow(panel)


#-------------------------------------------------------------------------------
# CONSORT flow diagram
#-------------------------------------------------------------------------------
# ---- Enrolment counts -------------------------------------------------------
# Real study flow: nobody was excluded after
# consenting, so there is no "Excluded" box between Consented and
# Randomised.
n_invited    <- 2160
n_consented  <- n_total
n_randomised <- n_consented

# Short description of what each arm actually received, per the trial protocol.
# Will be changed with real allocated arms and correct names after unblinding
arm_descriptions <- c(
  V1_control              = "Blinded (control)", # Current formulation (control)",
  V2_sentence             = "Blinded", # Added sentence about when it is okay to participate in activities or go to work",
  V3_definitions          = "Blinded", # Added definitions of key terms",
  V4_sentence_definitions = "Blinded" # Added sentence about activities/work + added definitions of key terms"
)

# Real allocation counts per arm (via the blinded file -> arm mapping in
# final_scripts/1_organize_data.R for the
# true per-file randomized totals). `panel` now contains one row per
# RANDOMISED participant per arm (1_organize_data.R pads the real completer
# rows back up to n_allocated_by_arm with placeholder "missing" rows, blank
# except participant_id/arm/included).
n_allocated_by_arm <- readRDS("panel_enrolment.rds")
n_analysed_by_arm  <- panel |>
  filter(!is.na(answer_time_ms)) |>
  count(arm) |>
  deframe() |>
  (\(x) x[names(arm_descriptions)])()

flow <- tibble(
  arm         = factor(names(arm_descriptions), levels = names(arm_descriptions)),
  n_allocated = as.integer(n_allocated_by_arm[names(arm_descriptions)]),
  n_analysed  = as.integer(n_analysed_by_arm)
) |>
  mutate(n_missing_outcome = n_allocated - n_analysed)
stopifnot(sum(flow$n_allocated) == n_randomised)

arm_short <- flow$arm |>
  as.character() |>
  str_extract("^V[0-9]+")

arm_desc <- str_wrap(arm_descriptions[as.character(flow$arm)], width = 16)

n_arms  <- nrow(flow)
x_col   <- seq(0, by = 7, length.out = n_arms)
x_mid   <- mean(x_col)
box_w   <- 2.6
box_h   <- 1.5
x_label <- min(x_col) - box_w / 2 - 2.2

# Enrolment segment: Invited -> Consented -> Randomised (nobody excluded).
y_invited    <- 15.5
y_consented  <- 13
y_randomised <- 10
y_alloc      <- 7
y_analysis   <- 3.5

box_invited <- tibble(
  xmin = x_mid - box_w, xmax = x_mid + box_w,
  ymin = y_invited - box_h / 2, ymax = y_invited + box_h / 2,
  label = paste0("Invited\n(N = ", n_invited, ")")
)

box_consented <- tibble(
  xmin = x_mid - box_w, xmax = x_mid + box_w,
  ymin = y_consented - box_h / 2, ymax = y_consented + box_h / 2,
  label = paste0("Consented\n(N = ", n_consented, ")")
)

box_randomised <- tibble(
  xmin = x_mid - box_w, xmax = x_mid + box_w,
  ymin = y_randomised - box_h / 2, ymax = y_randomised + box_h / 2,
  label = paste0("Randomised\n(N = ", n_randomised, ")")
)

box_alloc <- tibble(
  xmin = x_col - box_w / 2, xmax = x_col + box_w / 2,
  ymin = y_alloc - box_h / 2, ymax = y_alloc + box_h / 2,
  label = paste0(arm_short, "\n(n = ", flow$n_allocated, ")")
)

box_analysed <- tibble(
  xmin = x_col - box_w / 2, xmax = x_col + box_w / 2,
  ymin = y_analysis - box_h / 2, ymax = y_analysis + box_h / 2,
  label = paste0(arm_short, "\n(n = ", flow$n_analysed, ")")
)

# Midpoint of each arm's Allocation -> Analysis arrow, where a per-arm
# "Excluded" (missing outcome data) box branches off to the side.
y_mid_arm <- mean(c(y_alloc - box_h / 2, y_analysis + box_h / 2))
excl_arm_w <- 2.4

box_excl_arm <- tibble(
  xmin = x_col + box_w / 2 + 0.4, xmax = x_col + box_w / 2 + 0.4 + excl_arm_w,
  ymin = y_mid_arm - box_h / 2.2, ymax = y_mid_arm + box_h / 2.2,
  label = paste0("Missing data\n(n = ", flow$n_missing_outcome, ")")
)

boxes <- bind_rows(
  mutate(box_invited,    type = "main"),
  mutate(box_consented,  type = "main"),
  mutate(box_randomised, type = "main"),
  mutate(box_alloc,      type = "main"),
  mutate(box_analysed,   type = "main"),
  mutate(box_excl_arm,   type = "excl")
)

row_labels <- tibble(
  x = x_label,
  y = c(mean(c(y_invited, y_consented)), y_randomised, y_alloc, y_analysis),
  label = c("Enrolment", "Randomisation", "Allocation", "Analysis")
)

desc_labels <- tibble(
  x = x_col,
  y = y_analysis - box_h / 2 - 1.1,
  label = arm_desc
)

arrows_enrol <- tibble(
  x = x_mid,        y = c(y_invited - box_h / 2, y_consented - box_h / 2),
  xend = x_mid,      yend = c(y_consented + box_h / 2, y_randomised + box_h / 2)
)
arrows_split <- tibble(
  x = x_mid, y = y_randomised - box_h / 2,
  xend = x_col, yend = y_alloc + box_h / 2
)
arrows_down <- tibble(
  x = x_col, y = y_alloc - box_h / 2,
  xend = x_col, yend = y_analysis + box_h / 2
)
arrows_excl_arm <- tibble(
  x = x_col,                       y = y_mid_arm,
  xend = x_col + box_w / 2 + 0.4,  yend = y_mid_arm
)

consort_plot <- ggplot() +
  geom_segment(
    data = bind_rows(arrows_enrol, arrows_split, arrows_down),
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.4
  ) +
  geom_segment(
    data = arrows_excl_arm,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.4
  ) +
  geom_rect(
    data = boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, linetype = type),
    fill = "white", colour = "black"
  ) +
  geom_text(
    data = filter(boxes, type == "main"),
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 3, lineheight = 0.95
  ) +
  geom_text(
    data = filter(boxes, type == "excl"),
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 2.5, lineheight = 0.95
  ) +
  geom_text(
    data = row_labels,
    aes(x = x, y = y, label = label),
    fontface = "italic", size = 3.2, hjust = 0
  ) +
  geom_text(
    data = desc_labels,
    aes(x = x, y = y, label = label),
    fontface = "italic", size = 2.3, colour = "grey30",
    hjust = 0.5, vjust = 1, lineheight = 0.9
  ) +
  scale_linetype_manual(values = c(main = "solid", excl = "dashed"), guide = "none") +
  labs(
    title = "CONSORT flow chart",
    caption = str_wrap(
      "Missing data per arm: missing outcome data (did not answer the survey).",
      width = 90
    )
  ) +
  coord_cartesian(clip = "off", ylim = c(0.2, 16.5)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    plot.caption = element_text(hjust = 0.5, size = 9, margin = margin(t = 12))
  )

consort_plot

ggsave(
  file.path(results_dir, "consort_flowchart.png"),
  consort_plot, width = 9, height = 6.6, dpi = 150, bg = "white"
)

#-------------------------------------------------------------------------------
# Table 1: Baseline characteristics and behaviors
#-------------------------------------------------------------------------------
library(tidyr)
library(purrr)
library(forcats)
library(readr)

# Table 1 describes the whole randomized/allocated sample (not just the
# analysed subset), so participants with missing outcome data show up as an 
# explicit "Missing" category per characteristic below, rather than being 
# silently excluded.
baseline_pop <- panel |>
  rename("Age group" = "age_group",
         "Gender" = "gender",
         "Region" = "region",
         "Education" = "education",
         "Employed" = "employed",
         "Health literacy" = "health_literacy",
         "Can work from home" = "can_work_from_home",
         "Baseline behaviour" = "baseline_behaviour",
         "Contact high risk patients" = "contact_high_risk_patients")

baseline_vars <- c(
  "Age group", "Gender", "Region", "Education", "Health literacy", "Employed",
  "Contact high risk patients", "Can work from home", "Baseline behaviour"
)

# These three variables are skip-logic-routed (asked only of employed
# participants. For every other baseline variable, *everyone* is asked the 
# question, so any NA there is a dropout (missing_flag), not a routing outcome.
routed_vars <- c("Baseline behaviour", "Contact high risk patients", "Can work from home")

# Replace NA with an explicit category, distinguishing *why* it's missing:
# "Missing" for dropouts (no outcome data at all -- see simulate_panel_test.R),
# vs. "Not applicable" for structural skip logic (work-related questions asked
# only of the employed). "Not applicable" is only ever assigned for routed
# variables -- any other NA (which, given the simulation, shouldn't occur
# outside dropouts, but is guarded against here) falls back to "Missing"
# rather than being silently mislabelled as routing.
prep_var <- function(x, missing_flag, is_routed) {
  # Some baseline variables (e.g. age_group, education) are ordered factors;
  # dropping to plain character/factor here keeps the "level" column a
  # consistent type once map_dfr() stacks different variables' summaries.
  lvls <- levels(x)
  chr  <- as.character(x)
  chr  <- case_when(
    !is.na(chr)          ~ chr,
    missing_flag          ~ "Missing",
    is_routed             ~ "Not applicable",
    TRUE                  ~ "Missing"
  )
  factor(chr, levels = c(lvls, "Not applicable", "Missing"))
}

summarise_by_arm <- function(data, var) {
  data |>
    mutate(level = prep_var(.data[[var]], is.na(primary_correct_count), var %in% routed_vars)) |>
    count(arm, level, .drop = FALSE) |>
    group_by(arm) |>
    mutate(pct = 100 * n / sum(n)) |>
    ungroup() |>
    mutate(variable = var, cell = sprintf("%d (%.1f%%)", n, pct)) |>
    select(variable, level, arm, cell)
}

summarise_total <- function(data, var) {
  data |>
    mutate(level = prep_var(.data[[var]], is.na(primary_correct_count), var %in% routed_vars)) |>
    count(level, .drop = FALSE) |>
    mutate(pct = 100 * n / sum(n)) |>
    mutate(variable = var, Total = sprintf("%d (%.1f%%)", n, pct)) |>
    select(variable, level, Total)
}

by_arm    <- map_dfr(baseline_vars, ~ summarise_by_arm(baseline_pop, .x)) |>
  pivot_wider(names_from = arm, values_from = cell, values_fill = "0 (0.0%)")
total_col <- map_dfr(baseline_vars, ~ summarise_total(baseline_pop, .x))

table1 <- by_arm |>
  left_join(total_col, by = c("variable", "level")) |>
  arrange(match(variable, baseline_vars), level) |>
  rename(Characteristic = variable, Level = level) |>
  # "Not applicable" (skip logic) is only a real answer option for the
  # routed variables -- drop that row entirely for the others, where it's
  # structurally always 0.
  filter(Level != "Not applicable" | Characteristic %in% routed_vars)

# "Reasons for working from home" is a second layer of skip logic nested
# inside "Baseline behavior": it's only asked of participants who answered
# "Worked from home" there (see reasons_home_office_* in clean_panel_test.R).
# It's also a select-all-that-apply item (the four reasons aren't mutually
# exclusive -- a participant can endorse several), so it can't be summarized
# with prep_var()/summarize_by_arm() like the single-answer variables above:
# percentages are not expected to sum to 100% across its rows, and it's
# handled as its own block below instead, then appended to table1.
wfh_items <- tibble(
  label = c(
    "Too sick to work",
    "Concern about infecting others",
    "Concern that colleagues would dislike attendance",
    "Planned to work from home anyway"
  ),
  col = c(
    "reasons_home_office_too_sick",
    "reasons_home_office_infect",
    "reasons_home_office_stigma",
    "reasons_home_office_anyway"
  )
)
wfh_large_extent <- c("To a large extent", "To a very large extent")

# Eligibility is driven by "Baseline behavior" (not by employment status
# like the other routed_vars), so a dropout's NA there must not be read as
# "not eligible" -- treat dropouts (missing_flag) separately from genuine
# non-"Worked from home" answers.
wfh_eligible <- baseline_pop$`Baseline behaviour` == "Worked from home"
wfh_eligible[is.na(wfh_eligible)] <- FALSE
wfh_missing  <- is.na(baseline_pop$primary_correct_count)
wfh_arm_n    <- table(baseline_pop$arm)
wfh_total_n  <- nrow(baseline_pop)

# Among eligible (worked-from-home) participants, a handful did not answer
# a given reason item (genuine item non-response, distinct from both
# dropout and skip logic) -- these also count as "Missing" for that item,
# per-item, rather than being folded into the shared "Not applicable" row.
summarise_wfh_item <- function(data, col, label) {
  endorsed <- wfh_eligible & !wfh_missing & as.character(data[[col]]) %in% wfh_large_extent

  by_arm <- tibble(arm = data$arm, endorsed) |>
    group_by(arm) |>
    summarise(n = sum(endorsed), .groups = "drop") |>
    mutate(N = as.integer(wfh_arm_n[as.character(arm)]),
           cell = sprintf("%d (%.1f%%)", n, 100 * n / N)) |>
    select(arm, cell) |>
    pivot_wider(names_from = arm, values_from = cell)

  total_n <- sum(endorsed)
  by_arm |>
    mutate(Characteristic = "Reasons for working from home",
           Level = label,
           Total = sprintf("%d (%.1f%%)", total_n, 100 * total_n / wfh_total_n))
}

wfh_yes_rows <- map2_dfr(wfh_items$col, wfh_items$label,
                          ~ summarise_wfh_item(baseline_pop, .x, .y))

# Shared "Not applicable" (didn't work from home) / "Missing" (dropout)
# rows, computed once since eligibility and dropout status don't vary
# across the four reason items.
wfh_cat <- case_when(
  wfh_missing    ~ "Missing",
  !wfh_eligible  ~ "Not applicable",
  TRUE           ~ NA_character_
)

wfh_na_rows <- tibble(arm = baseline_pop$arm, cat = wfh_cat) |>
  filter(!is.na(cat)) |>
  count(arm, cat, name = "n") |>
  complete(arm = levels(baseline_pop$arm), cat = c("Not applicable", "Missing"),
           fill = list(n = 0)) |>
  mutate(N = as.integer(wfh_arm_n[as.character(arm)]),
         cell = sprintf("%d (%.1f%%)", n, 100 * n / N)) |>
  select(arm, cat, cell) |>
  pivot_wider(names_from = arm, values_from = cell) |>
  rename(Level = cat)

wfh_na_total <- tibble(cat = wfh_cat) |>
  filter(!is.na(cat)) |>
  count(cat, name = "n") |>
  complete(cat = c("Not applicable", "Missing"), fill = list(n = 0)) |>
  mutate(Total = sprintf("%d (%.1f%%)", n, 100 * n / wfh_total_n)) |>
  select(cat, Total) |>
  rename(Level = cat)

wfh_na_rows <- wfh_na_rows |>
  left_join(wfh_na_total, by = "Level") |>
  mutate(Characteristic = "Reasons for working from home")

wfh_rows <- bind_rows(wfh_yes_rows, wfh_na_rows) |>
  mutate(Level = factor(Level, levels = c(wfh_items$label, "Not applicable", "Missing"))) |>
  arrange(Level) |>
  mutate(Level = as.character(Level)) |>
  select(Characteristic, Level, all_of(names(table1)[-(1:2)]))

# Placed right after "Baseline behavior", the question it's nested under.
table1 <- bind_rows(table1, wfh_rows)
baseline_vars <- c(baseline_vars, "Reasons for working from home")
routed_vars   <- c(routed_vars, "Reasons for working from home")

table1

write_csv(table1, file.path(results_dir, "table1_baseline_characteristics.csv"))

# ---- Publication-ready version (gt) ----
# Saved as .docx/.rtf so it can be pasted directly into a manuscript, plus
# .html for quick viewing.
library(gt)

col_n   <- setNames(flow$n_allocated, as.character(flow$arm))
total_n <- sum(flow$n_allocated)

table1_gt <- table1 |>
  gt(groupname_col = "Characteristic", rowname_col = "Level") |>
  tab_header(
    title    = "Table 1. Baseline characteristics and behaviours",
    subtitle = "By study arm, among randomised/allocated participants"
  ) |>
  cols_label(
    V1_control              = html(sprintf("V1<br>(n = %d)", col_n["V1_control"])),
    V2_sentence              = html(sprintf("V2<br>(n = %d)", col_n["V2_sentence"])),
    V3_definitions           = html(sprintf("V3<br>(n = %d)", col_n["V3_definitions"])),
    V4_sentence_definitions  = html(sprintf("V4<br>(n = %d)", col_n["V4_sentence_definitions"])),
    Total                    = html(sprintf("Total<br>(n = %d)", total_n))
  ) |>
  tab_spanner(
    label   = "Study arm, n (%)",
    columns = c(V1_control, V2_sentence, V3_definitions, V4_sentence_definitions)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_text(align = "center"),
    locations = cells_body(columns = c(V1_control, V2_sentence, V3_definitions, V4_sentence_definitions, Total))
  ) |>
  tab_style(
    style     = cell_text(align = "center"),
    locations = cells_column_labels(columns = c(V1_control, V2_sentence, V3_definitions, V4_sentence_definitions, Total))
  ) |>
  cols_align(align = "left", columns = Level) |>
  tab_footnote(
    footnote = "V1 = control (current formulation); V2 = added sentence about when it is okay to participate in activities or go to work; V3 = added definitions of key terms; V4 = added sentence and definitions."
  ) |>
  tab_footnote(
    footnote = "\"Missing\" = no outcome/survey data at all (dropped out before completing the primary-outcome questions; see the CONSORT flow chart's per-arm \"Missing data\" counts)."
  ) |>
  tab_footnote(
    footnote = "\"Not applicable\" reflects skip-logic among participants who did respond (e.g. work-related items asked only of employed participants, or behaviour-specific reasons asked only of participants reporting that behaviour)."
  ) |>
  tab_footnote(
    footnote  = "Participants could endorse more than one reason, so percentages need not sum to 100%; \"Missing\" here also includes eligible (worked-from-home) participants who did not answer that specific reason item, in addition to dropouts.",
    locations = cells_row_groups(groups = "Reasons for working from home")
  ) |>
  tab_options(
    table.font.size             = px(12),
    heading.title.font.size     = px(14),
    heading.subtitle.font.size  = px(12),
    column_labels.font.weight   = "bold",
    table.border.top.style      = "solid",
    table.border.bottom.style   = "solid"
  ) |>
  opt_table_font(font = "Times New Roman")

table1_gt

gtsave(table1_gt, file.path(results_dir, "table1_baseline_characteristics.docx"))
gtsave(table1_gt, file.path(results_dir, "table1_baseline_characteristics.rtf"))
gtsave(table1_gt, file.path(results_dir, "table1_baseline_characteristics.html"))

#-------------------------------------------------------------------------------


