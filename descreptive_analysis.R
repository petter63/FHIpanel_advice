#-------------------------------------------------------------------------------
# Prosjekt: FHI-panel effekt av smittevernraad
#-------------------------------------------------------------------------------
# Beskrivelse av test-data

panel_test <- readRDS("panel_test.rds")
str(panel_test)

#-------------------------------------------------------------------------------
# CONSORT flytdiagram
#-------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)

# Each run gets its own timestamped subfolder under results/, so output from
# different runs is kept separate and it's clear when a given result was produced.
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

n_total <- nrow(panel_test)

# Number randomised/excluded/analysed per arm.
# Exclusion here: failed the attention check (attention_check_pass == FALSE).
flow <- panel_test |>
  count(arm, attention_check_pass) |>
  tidyr::pivot_wider(names_from = attention_check_pass, values_from = n, values_fill = 0) |>
  rename(n_excluded = `FALSE`, n_analysed = `TRUE`) |>
  mutate(n_allocated = n_excluded + n_analysed) |>
  arrange(arm)

arm_short <- flow$arm |>
  as.character() |>
  str_extract("^V[0-9]+")

# Short description of what each arm actually received, per the trial protocol.
arm_descriptions <- c(
  V1_control              = "Current formulation (control)",
  V2_sentence             = "Added sentence about when it is okay to participate in activities or go to work",
  V3_definitions          = "Added definitions of key terms",
  V4_sentence_definitions = "Added sentence about activities/work + added definitions of key terms"
)
arm_desc <- str_wrap(arm_descriptions[as.character(flow$arm)], width = 16)

n_arms  <- nrow(flow)
x_col   <- seq(0, by = 7, length.out = n_arms)
box_w   <- 2.6
box_h   <- 1.5
excl_w  <- 2.2
x_label <- min(x_col) - box_w / 2 - 2.2

box_randomised <- tibble(
  xmin = mean(x_col) - box_w, xmax = mean(x_col) + box_w,
  ymin = 10 - box_h / 2,      ymax = 10 + box_h / 2,
  label = paste0("Randomised\n(N = ", n_total, ")")
)

box_alloc <- tibble(
  xmin = x_col - box_w / 2, xmax = x_col + box_w / 2,
  ymin = 7 - box_h / 2,     ymax = 7 + box_h / 2,
  label = paste0(arm_short, "\n(n = ", flow$n_allocated, ")")
)

box_excl <- tibble(
  xmin = x_col + box_w / 2 + 0.4, xmax = x_col + box_w / 2 + 0.4 + excl_w,
  ymin = 7 - box_h / 2,           ymax = 7 + box_h / 2,
  label = paste0("Excl.\n(n = ", flow$n_excluded, ")")
)

box_analysed <- tibble(
  xmin = x_col - box_w / 2, xmax = x_col + box_w / 2,
  ymin = 3.5 - box_h / 2,   ymax = 3.5 + box_h / 2,
  label = paste0(arm_short, "\n(n = ", flow$n_analysed, ")")
)

boxes <- bind_rows(
  mutate(box_randomised, type = "main"),
  mutate(box_alloc,      type = "main"),
  mutate(box_excl,       type = "excl"),
  mutate(box_analysed,   type = "main")
)

row_labels <- tibble(
  x = x_label,
  y = c(10, 7, 3.5),
  label = c("Randomisation", "Allocation", "Analysis")
)

desc_labels <- tibble(
  x = x_col,
  y = 3.5 - box_h / 2 - 1.1,
  label = arm_desc
)

arrows_split <- tibble(
  x = mean(x_col), y = 10 - box_h / 2,
  xend = x_col,     yend = 7 + box_h / 2
)
arrows_excl <- tibble(
  x = x_col + box_w / 2,          y = 7,
  xend = x_col + box_w / 2 + 0.4, yend = 7
)
arrows_down <- tibble(
  x = x_col, y = 7 - box_h / 2,
  xend = x_col, yend = 3.5 + box_h / 2
)

consort_plot <- ggplot() +
  geom_segment(
    data = bind_rows(arrows_split, arrows_down),
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.4
  ) +
  geom_segment(
    data = arrows_excl,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.4
  ) +
  geom_rect(
    data = boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, linetype = type),
    fill = "white", colour = "black"
  ) +
  geom_text(
    data = boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 3, lineheight = 0.95
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
    caption = "Exclusion: failed the attention check"
  ) +
  coord_cartesian(clip = "off", ylim = c(0.2, 11)) +
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
# Table 2: Baseline characteristics and behaviours
#-------------------------------------------------------------------------------
library(tidyr)
library(purrr)
library(forcats)
library(readr)

# Restrict to the analysed sample (excludes those who failed the attention
# check), consistent with the "Analysed" count in the CONSORT flow chart.
analysed <- panel_test |> filter(attention_check_pass)

baseline_vars <- c(
  "age_group", "gender", "region", "education", "employed", "health_literacy",
  "baseline_behaviour", "contact_high_risk_patients", "can_work_from_home"
)

# Replace NA with an explicit "Not applicable" category (e.g. the work-related
# questions are only asked of those who are employed), so denominators stay clear.
prep_var <- function(x) {
  if (is.factor(x)) {
    fct_na_value_to_level(x, "Not applicable")
  } else {
    replace_na(x, "Not applicable")
  }
}

summarise_by_arm <- function(data, var) {
  data |>
    mutate(level = prep_var(.data[[var]])) |>
    count(arm, level) |>
    group_by(arm) |>
    mutate(pct = 100 * n / sum(n)) |>
    ungroup() |>
    mutate(variable = var, cell = sprintf("%d (%.1f%%)", n, pct)) |>
    select(variable, level, arm, cell)
}

summarise_total <- function(data, var) {
  data |>
    mutate(level = prep_var(.data[[var]])) |>
    count(level) |>
    mutate(pct = 100 * n / sum(n)) |>
    mutate(variable = var, Total = sprintf("%d (%.1f%%)", n, pct)) |>
    select(variable, level, Total)
}

by_arm    <- map_dfr(baseline_vars, ~ summarise_by_arm(analysed, .x)) |>
  pivot_wider(names_from = arm, values_from = cell, values_fill = "0 (0.0%)")
total_col <- map_dfr(baseline_vars, ~ summarise_total(analysed, .x))

table2 <- by_arm |>
  left_join(total_col, by = c("variable", "level")) |>
  arrange(match(variable, baseline_vars), level) |>
  rename(Characteristic = variable, Level = level)

table2

write_csv(table2, file.path(results_dir, "table2_baseline_characteristics.csv"))

# ---- Publication-ready version (gt) ----
# Saved as .docx/.rtf so it can be pasted directly into a manuscript, plus
# .html for quick viewing.
library(gt)

col_n   <- setNames(flow$n_analysed, as.character(flow$arm))
total_n <- sum(flow$n_analysed)

table2_gt <- table2 |>
  gt(groupname_col = "Characteristic", rowname_col = "Level") |>
  tab_header(
    title    = "Table 2. Baseline characteristics and behaviours",
    subtitle = "By study arm, among analysed participants"
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
    footnote = "Restricted to participants who passed the attention check."
  ) |>
  tab_footnote(
    footnote = "\"Not applicable\" reflects skip-logic (e.g. work-related items asked only of employed participants, or behaviour-specific reasons asked only of participants reporting that behaviour)."
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

table2_gt

gtsave(table2_gt, file.path(results_dir, "table2_baseline_characteristics.docx"))
gtsave(table2_gt, file.path(results_dir, "table2_baseline_characteristics.rtf"))
gtsave(table2_gt, file.path(results_dir, "table2_baseline_characteristics.html"))

#-------------------------------------------------------------------------------
