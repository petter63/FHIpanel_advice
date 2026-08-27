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

# Short description of what each arm actually received.
arm_descriptions <- c(
  V1_explanations          = "Advice with explanations",
  V2_explanations_sentence = "Advice with explanations + one-sentence summary",
  V3_sentence_only         = "Advice with one-sentence summary only",
  V4_control               = "Standard advice (control)"
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
    title = "CONSORT flow diagram",
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

