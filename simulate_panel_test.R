## Simulate a single test data frame with exactly the same variables as the
## four real result files in data/ ("Smittevern re[s]sult[s]_{W,X,Y,Z}.xlsx"),
## plus one added "arm" variable identifying which of the four files/arms
## each simulated respondent belongs to.
##
## Each real file is one arm's raw, uncleaned export: metadata columns
## prefixed with "$", categorical answers coded as numeric strings,
## multi-select questions stored as separate 0/1 checkbox columns, and
## survey skip logic producing NA where a question was not shown. The raw
## files themselves carry no arm label -- which file a respondent is in
## *is* their arm -- so that information is added here as an explicit
## "arm" column once the four arms are combined into one data frame.
##
## Column names/values below (probabilities, coding, skip logic) were read
## directly off the four real files, using each file's canonical column
## names (the real files have a couple of minor per-file naming
## inconsistencies -- e.g. a typo'd column and inconsistent casing on
## Scenario_3 items in "Smittevern results_Y.xlsx" -- which are normalised
## away here since this script targets one consistent variable set).
##
## The scoring rules for the comprehension items and the two checkbox
## outcomes are taken from "data/Smittevern_RCT_codebook and data
## dictionary.xlsx" (sheet "Detailed variable key" / "Derived outcomes"):
##  - Scenario_*.Careful is *always* correct at the high end (4-5), and
##    Scenario_*.Normal is *always* correct at the low end (1-2),
##    regardless of scenario -- only the Delay/Home item's correct
##    direction flips by scenario (see `scenario_spec$correct_high`).
##  - High_risk checkbox: correct = options {1, 2, 4, 6} selected and
##    {3, 5} not selected (option suffixes follow the questionnaire's
##    option numbers, not worksheet order -- .6 is the healthy 2-month-old
##    and .5 is "none of the above").

rm(list = ls())

library(tidyverse)

set.seed(6274)

## ---- Arm / file setup ---------------------------------------------------
## Assumed mapping to the trial protocol (adjust if the true mapping
## differs -- the raw files carry no explicit arm label). The codebook's
## own "Review before analysis" sheet flags this as a CRITICAL open issue
## ("Randomized arm absent" -> "Merge verified allocation data"), so this
## mapping remains a placeholder until the real randomisation log is
## merged in:
##   W -> V1 (control), X -> V2 (+ sentence), Y -> V3 (+ definitions),
##   Z -> V4 (+ sentence & definitions)
arm_files <- tibble(
  file_label       = c("W", "X", "Y", "Z"),
  arm              = c("V1_control", "V2_sentence", "V3_definitions", "V4_sentence_definitions"),
  n                = c(224, 218, 220, 214),
  p_correct        = c(0.55, 0.68, 0.65, 0.72),
  has_explanations = c(FALSE, FALSE, TRUE, TRUE)
)

## ---- Helpers -------------------------------------------------------------

# Directional 5-point Likert item: with probability p_correct, the answer
# leans toward the "correct" end (1-2 or 4-5); otherwise uniform over 1-5.
sim_likert_directional <- function(n, correct_high, p_correct) {
  is_correct <- runif(n) < p_correct
  vapply(seq_len(n), function(i) {
    if (is_correct[i]) {
      if (correct_high) sample(4:5, 1) else sample(1:2, 1)
    } else {
      sample(1:5, 1)
    }
  }, integer(1))
}

# Sample a coded categorical item (1..length(probs)) and return it as
# character, with optional sparse missingness (mirrors the small amount of
# genuine item non-response present in the real exports).
sim_coded <- function(n, probs, p_na = 0) {
  codes <- as.character(sample(seq_along(probs), n, replace = TRUE, prob = probs))
  if (p_na > 0) codes[runif(n) < p_na] <- NA_character_
  codes
}

# Multi-select checkbox block (0/1 per option, stored as character to match
# the raw export).
sim_checkbox_block <- function(n, cols, true_opts, p_true_hit, p_false_hit) {
  false_opts <- setdiff(seq_along(cols), true_opts)
  mat <- matrix(0L, nrow = n, ncol = length(cols), dimnames = list(NULL, cols))
  for (i in seq_len(n)) {
    mat[i, true_opts]  <- as.integer(runif(length(true_opts))  < p_true_hit[i])
    mat[i, false_opts] <- as.integer(runif(length(false_opts)) < p_false_hit[i])
  }
  as_tibble(mat) |> mutate(across(everything(), as.character))
}

# The 5 comprehension scenarios: each has 3 directional items + an
# "enough info" Likert item. Column name stems and casing follow each real
# file's canonical (most common) naming, which is inconsistent across
# scenarios -- e.g. "Scenario_1.careful" (lowercase) vs "Scenario_4.Careful"
# (uppercase).
scenario_spec <- tibble(
  scenario     = 1:5,
  item1_name   = c("Delay", "Delay", "Home", "Delay", "Home"),
  correct_high = c(TRUE, TRUE, FALSE, FALSE, FALSE),
  item23_case  = c("lower", "lower", "upper", "upper", "upper")
)

sim_scenario_block <- function(n, p_correct) {
  map_dfc(seq_len(nrow(scenario_spec)), function(s) {
    spec <- scenario_spec[s, ]
    # Per the codebook: Careful is always correct-high, Normal is always
    # correct-low, independent of scenario; only Delay/Home (item1) flips.
    item1 <- sim_likert_directional(n, spec$correct_high, p_correct)
    item2 <- sim_likert_directional(n, TRUE,  p_correct * 0.9)
    item3 <- sim_likert_directional(n, FALSE, p_correct * 0.9)
    enough_info <- sample(1:5, n, replace = TRUE, prob = c(0.05, 0.08, 0.12, 0.35, 0.40))

    use_lower <- spec$item23_case == "lower"
    careful_name <- if (use_lower) "careful" else "Careful"
    normal_name  <- if (use_lower) "normal"  else "Normal"

    tibble(
      !!paste0("Scenario_", spec$scenario, ".", spec$item1_name) := as.character(item1),
      !!paste0("Scenario_", spec$scenario, ".", careful_name)    := as.character(item2),
      !!paste0("Scenario_", spec$scenario, ".", normal_name)     := as.character(item3),
      !!paste0("Scenario_", spec$scenario, "_enough_info")       := as.character(enough_info)
    )
  })
}

## ---- Per-arm simulation ---------------------------------------------------

simulate_one_arm <- function(arm, n, p_correct, has_explanations,
                              submission_ids, created) {

  ## Employment: single-select checkbox block (Employment.1/2/3)
  employment_choice <- sample(1:3, n, replace = TRUE, prob = c(0.805, 0.167, 0.029))
  employed <- employment_choice == 1

  ## Work-related items, asked only of the employed
  close_contact  <- ifelse(employed, sim_coded(n, c(0.237, 0.685, 0.078)), NA)
  home_office    <- ifelse(employed, sim_coded(n, c(0.554, 0.438, 0.009)), NA)
  last_time_sick <- ifelse(employed,
                            sim_coded(n, c(0.152, 0.388, 0.229, 0.036, 0.173, 0.021, 0.001)),
                            NA)

  is_wfh_episode       <- !is.na(last_time_sick) & last_time_sick == "3"
  is_sick_note_episode <- !is.na(last_time_sick) & last_time_sick == "5"

  last_time_home_office <- ifelse(is_wfh_episode, sim_coded(n, c(0.944, 0.056)), NA)

  reason_cols_wfh <- c("Reasons_home_office.too_sick", "Reasons_home_office.infect",
                        "Reasons_home_office.stigma", "Reasons_home_office.home_office_anyway")
  reasons_home_office <- map_dfc(reason_cols_wfh, function(col) {
    asked <- is_wfh_episode & runif(n) > 0.08  # ~8% real-world non-response
    tibble(!!col := as.character(ifelse(asked, sample(1:5, n, replace = TRUE), NA)))
  })

  sick_note     <- ifelse(is_sick_note_episode, sim_coded(n, c(0.95, 0.05)), NA)
  has_sick_note <- !is.na(sick_note) & sick_note == "1"

  reason_cols_sick <- c("Reasons_sick_note.too_sick", "Reasons_sick_note.infect",
                         "Reasons_sick_note.stigma")
  reasons_sick_note <- map_dfc(reason_cols_sick, function(col) {
    tibble(!!col := as.character(ifelse(has_sick_note, sample(1:5, n, replace = TRUE), NA)))
  })

  ## Comprehension scenarios (primary outcome) - the only block that varies
  ## systematically by arm, via p_correct.
  scenario_items <- sim_scenario_block(n, p_correct)

  ## Secondary outcomes: some vary by has_explanations, as in the trial's
  ## hypothesis that added definitions/sentences improve comprehension.
  p_goal_correct <- if (has_explanations) c(0.02, 0.82, 0.03, 0.13) else c(0.02, 0.70, 0.04, 0.24)

  p_true_hit_infant  <- if (has_explanations) 0.80 else 0.55
  p_false_hit_infant <- if (has_explanations) 0.12 else 0.30
  infant_selection <- sim_checkbox_block(
    n, paste0("Infant.", 1:5), true_opts = c(1, 2, 4),
    p_true_hit = rep(p_true_hit_infant, n), p_false_hit = rep(p_false_hit_infant, n)
  )

  p_true_hit_risk  <- if (has_explanations) 0.78 else 0.50
  p_false_hit_risk <- if (has_explanations) 0.15 else 0.35
  risk_selection <- sim_checkbox_block(
    # Correct = {1, 2, 4, 6} selected; {3, 5} not selected (see codebook note
    # above on option suffixes vs. worksheet order).
    n, paste0("High_risk.", 1:6), true_opts = c(1, 2, 4, 6),
    p_true_hit = rep(p_true_hit_risk, n), p_false_hit = rep(p_false_hit_risk, n)
  )

  ## ---- Assemble: "arm" first, then the real files' variables in their
  ## original column order ----
  out <- tibble(
    arm               = arm,
    `$submission_id` = submission_ids,
    `$created`        = created,
    Gender            = sim_coded(n, c(0.570, 0.424, 0.002, 0.003)),
    Age               = sim_coded(n, c(0.024, 0.141, 0.182, 0.257, 0.355, 0.041)),
    Municipality      = sim_coded(n, c(0.071, 0.121, 0.165, 0.064, 0.063, 0.047, 0.030,
                                        0.043, 0.087, 0.092, 0.045, 0.086, 0.049, 0.026, 0.011)),
    Education         = sim_coded(n, c(0.023, 0.191, 0.123, 0.321, 0.342)),
    research_participation_experience = sim_coded(n, c(0.510, 0.156, 0.334), p_na = 0.05),
    Employment.1      = as.character(as.integer(employment_choice == 1)),
    Employment.2      = as.character(as.integer(employment_choice == 2)),
    Employment.3      = as.character(as.integer(employment_choice == 3)),
    Close_contact_patients = close_contact,
    Home_office        = home_office,
    Last_time_sick     = last_time_sick,
    Last_time_home_office = last_time_home_office
  ) |>
    bind_cols(reasons_home_office) |>
    bind_cols(tibble(Sick_note = sick_note)) |>
    bind_cols(reasons_sick_note) |>
    bind_cols(tibble(Health_literacy = sim_coded(n, c(0.022, 0.147, 0.638, 0.193)))) |>
    bind_cols(scenario_items) |>
    bind_cols(tibble(
      Main_aim_advice = sim_coded(n, p_goal_correct),
      Confident_understanding = sim_coded(n, c(0.420, 0.553, 0.019, 0.008))
    )) |>
    bind_cols(infant_selection) |>
    bind_cols(risk_selection) |>
    bind_cols(tibble(
      Trustworthy = sim_coded(n, c(0.038, 0.013, 0.056, 0.489, 0.405)),
      Seek_info   = sim_coded(n, c(0.450, 0.505, 0.030, 0.015)),
      Share_info  = sim_coded(n, c(0.445, 0.489, 0.037, 0.030)),
      `$answer_time_ms`     = pmin(7.2e7, round(rlnorm(n, meanlog = log(4.3e5), sdlog = 0.9))),
      `$forwarded_to_form`  = 608092
    ))

  out
}

## ---- Generate metadata (submission ids / timestamps) shared across arms --

n_total <- sum(arm_files$n)
date_start <- as.POSIXct("2026-04-27", tz = "UTC")
date_end   <- as.POSIXct("2026-05-30", tz = "UTC")

id_pool      <- sort(sample(45480000:46100000, n_total))
created_pool <- sort(date_start + runif(n_total) * as.numeric(date_end - date_start, units = "secs"))

## ---- Simulate each arm and combine into one data frame ----

arm_levels <- arm_files$arm

panel_test_raw <- list()
idx <- 1
for (i in seq_len(nrow(arm_files))) {
  n_i   <- arm_files$n[i]
  rows  <- idx:(idx + n_i - 1)
  panel_test_raw[[i]] <- simulate_one_arm(
    arm              = factor(arm_files$arm[i], levels = arm_levels),
    n                = n_i,
    p_correct        = arm_files$p_correct[i],
    has_explanations = arm_files$has_explanations[i],
    submission_ids   = id_pool[rows],
    created          = created_pool[rows]
  )
  idx <- idx + n_i
}

# panel_test_raw is a single data frame stacking all four simulated arms: it
# has exactly the same variables as the real "Smittevern *.xlsx" files
# (same names, coding, and skip logic), plus one added "arm" column
# identifying which of the four files/arms each row belongs to.
panel_test_raw <- bind_rows(panel_test_raw)
saveRDS(panel_test_raw, "panel_test_raw.rds")
