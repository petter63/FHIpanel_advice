## Read the four REAL, blinded raw survey-export files (one per study arm,
## "data/Smittevern re[s]sult[s]_{W,X,Y,Z}.xlsx") and clean them into the
## analysis-ready tibble structure the rest of the pipeline expects.
##
## Value labels and derived outcomes below come from
## "data/Smittevern_RCT_codebook and data dictionary.xlsx" (sheets
## "Detailed variable key" and "Derived outcomes"). A few things to flag:
##
##  - CRITICAL (per the codebook's own "Review before analysis" sheet):
##    the real, verified randomised allocation (which file/letter -> which
##    of V1-V4) is *not* in the raw export, and is not known to whoever
##    exported these files either -- the four files must be read and
##    analysed BLINDED to arm. `arm` below is therefore assigned completely
##    AT RANDOM, one arm per file, purely as a placeholder so the rest of
##    the pipeline (which expects an `arm` column) runs end to end. This
##    is NOT the real allocation and every arm-stratified result produced
##    from this file is meaningless until the real randomisation log is
##    obtained and merged in -- replace the random assignment below with
##    it at that point.
##  - The four files have a few per-file naming inconsistencies (a typo'd
##    column and inconsistent casing on the Scenario_3 items in
##    "Smittevern results_Y.xlsx"), normalised away below so the four
##    files share one consistent column set before being combined.
##  - None of the four files contains an "Included" column, nor any rows
##    for people who started but never finished the survey -- these
##    exports apparently only contain completed submissions. "included" is
##    therefore set to TRUE for everyone below (same placeholder treatment
##    as attention_check_pass), and there is no "missing outcome data"
##    (dropout) row in this data at all. This means `n_allocated` per arm
##    (needed for the CONSORT flow chart in descreptive_analysis.R) is NOT
##    recoverable from these files and must come from the real
##    randomisation/enrolment log instead.
##  - Health_literacy is coded 1 = very difficult ... 4 = very easy (the
##    codebook flags this direction as differing from an earlier protocol
##    description -- use the administered-form coding, as done here).
##  - High_risk correct options are {1, 2, 4, 6}; {3, 5} are incorrect
##    (suffixes follow the questionnaire's option numbers, not worksheet
##    order).
##  - Scenario_*_enough_info's value labels are NOT in a monotonic
##    agreement order in the raw export (codes are
##    1=Disagree, 2=Strongly disagree, 3=Strongly agree, 4=Agree,
##    5=Neither agree nor disagree) -- decoded below with the levels
##    reordered into their actual agreement order.

rm(list = ls())

library(tidyverse)
library(readxl)

## ---- 1. Read the four raw files and tag each with a BLINDED, random arm --

data_dir <- file.path("data")

arm_levels <- c("V1_control", "V2_sentence", "V3_definitions", "V4_sentence_definitions")

raw_files <- c(
  W = file.path(data_dir, "Smittevern reslults_W.xlsx"),
  X = file.path(data_dir, "Smittevern results_X.xlsx"),
  Y = file.path(data_dir, "Smittevern results_Y.xlsx"),
  Z = file.path(data_dir, "Smittevern results_Z.xlsx")
)

# Blinded placeholder allocation: which file gets which arm label is random,
# not the real (unknown, to us) randomization -- see header note. Seeded
# only so re-running this script doesn't reshuffle the placeholder mapping
# underneath an already-started analysis; it carries no other meaning.
set.seed(8421)
file_arm_map <- setNames(sample(arm_levels), names(raw_files))
file_arm_map

# True total randomised per FILE (not per arm -- the real randomisation log
# itself is still not available, see header note), supplied 2026-10-05:
# 260 (Z), 250 (Y), 251 (X), 255 (W). Mapped through the same blinded
# file -> arm assignment above so the CONSORT counts stay consistent with
# whatever that random mapping happens to be. If file_arm_map above ever
# changes (different seed, different sample() call, etc.) this stays
# correct automatically since it's keyed off the same object.
n_allocated_by_file <- c(W = 255, X = 251, Y = 250, Z = 260)
stopifnot(sum(n_allocated_by_file) == 1016)  # matches total consented/randomised
n_allocated_by_arm  <- setNames(n_allocated_by_file[names(file_arm_map)], file_arm_map)
n_allocated_by_arm  <- n_allocated_by_arm[arm_levels]   # V1-V4 order
n_allocated_by_arm

raw_list <- map(raw_files, ~ read_excel(.x, sheet = 1))

# Normalise the per-file naming inconsistencies (see header note) to one
# consistent column set before combining.
raw_list$Y <- raw_list$Y |>
  rename(
    research_participation_experience = research_participant_experience,
    Scenario_3.Careful                = Scenario_3.careful,
    Scenario_3.Normal                 = Scenario_3.normal
  )

# No file has an "Included" column (see header note) -- add the same
# placeholder used for attention_check_pass below (everyone treated as
# having met the inclusion criteria, since there's no raw-data equivalent
# to check against).
raw_list <- map(raw_list, ~ mutate(.x, Included = "1"))

panel_data <- imap(raw_list, function(df, file_label) {
  mutate(df, arm = file_arm_map[[file_label]])
}) |>
  bind_rows()

# One arm per file: there should be exactly as many distinct arms as files,
# and each arm should appear in exactly one file's worth of rows.
stopifnot(n_distinct(panel_data$arm) == length(raw_files))

panel_data |> count(arm)

## ---- Label decoding helpers -----------------------------------------------

# `labels` is a named character vector, e.g. c("1" = "Female", "2" = "Male").
# Decodes raw character codes to a factor with those labels, in label order
# (or in `level_order` order, if a different -- e.g. semantic -- order than
# the raw codes is needed, as for Scenario_*_enough_info below).
decode <- function(x, labels, ordered = FALSE, level_order = names(labels)) {
  factor(unname(labels[x]), levels = unname(labels[level_order]), ordered = ordered)
}

gender_labels     <- c("1" = "Female", "2" = "Male", "3" = "Other", "4" = "Prefer not to say")
age_labels        <- c("1" = "16-24", "2" = "25-34", "3" = "35-44", "4" = "45-54",
                        "5" = "55-64", "6" = "65 or older")
region_labels     <- c("1" = "\u00d8stfold", "2" = "Akershus", "3" = "Oslo", "4" = "Innlandet",
                        "5" = "Buskerud", "6" = "Vestfold", "7" = "Telemark", "8" = "Agder",
                        "9" = "Rogaland", "10" = "Vestland", "11" = "M\u00f8re og Romsdal",
                        "12" = "Tr\u00f8ndelag", "13" = "Nordland", "14" = "Troms", "15" = "Finnmark")
education_labels  <- c("1" = "Primary school or lower", "2" = "Upper secondary school",
                        "3" = "Vocational college", "4" = "University/college, 4 years or less",
                        "5" = "University/college, more than 4 years")
research_exp_labels <- c("1" = "Yes", "2" = "No", "3" = "Don't remember")
yes_no_dontknow_labels <- c("1" = "Yes", "2" = "No", "3" = "Don't know")
last_time_sick_labels <- c(
  "1" = "Not applicable / no respiratory symptoms in past 12 months / happened during a holiday or weekend",
  "2" = "Was physically at work",
  "3" = "Worked from home",
  "4" = "Got a doctor's sick note",
  "5" = "Used self-certified sick leave",
  "6" = "Was not working at that time",
  "7" = "Don't remember"
)
yes_no_labels     <- c("1" = "Yes", "2" = "No")
extent_labels     <- c("1" = "Not at all", "2" = "To a small extent", "3" = "To some extent",
                        "4" = "To a large extent", "5" = "To a very large extent")
health_lit_labels <- c("1" = "Very difficult", "2" = "Quite difficult",
                        "3" = "Quite easy", "4" = "Very easy")
likelihood_labels <- c("1" = "Very unlikely", "2" = "Unlikely", "3" = "Don't know",
                        "4" = "Likely", "5" = "Very likely")
main_aim_labels   <- c(
  "1" = "Protect yourself",
  "2" = "Protect people at higher risk of becoming severely ill",
  "3" = "Protect others around me",
  "4" = "Prevent spread of infection in the community as much as possible"
)
confidence_labels <- c("1" = "Very confident", "2" = "Fairly confident",
                        "3" = "Fairly unconfident", "4" = "Very unconfident")
trust_labels      <- c("1" = "Very little", "2" = "A little", "3" = "Neither little nor much",
                        "4" = "Much", "5" = "Very much")
agreement_labels  <- c("1" = "Strongly agree", "2" = "Agree", "3" = "Disagree", "4" = "Strongly disagree")

# Scenario_*_enough_info: raw codes are NOT in agreement order (see header
# note); reorder into Strongly disagree < Disagree < Neither < Agree <
# Strongly agree.
enough_info_labels <- c("1" = "Disagree", "2" = "Strongly disagree", "3" = "Strongly agree",
                         "4" = "Agree", "5" = "Neither agree nor disagree")
enough_info_order   <- c("2", "1", "5", "4", "3")

## ---- Comprehension scenario items -----------------------------------------

# Column names/casing are fixed per scenario (now normalised across files,
# see above); all three items per scenario share the same "likelihood"
# response scale.
scenario_cols <- tribble(
  ~scenario, ~item1,             ~item2,               ~item3,
  1,         "Scenario_1.Delay", "Scenario_1.careful", "Scenario_1.normal",
  2,         "Scenario_2.Delay", "Scenario_2.careful", "Scenario_2.normal",
  3,         "Scenario_3.Home",  "Scenario_3.Careful",  "Scenario_3.Normal",
  4,         "Scenario_4.Delay", "Scenario_4.Careful",  "Scenario_4.Normal",
  5,         "Scenario_5.Home",  "Scenario_5.Careful",  "Scenario_5.Normal"
)
# Correct direction for item1 (Delay/Home) only; item2 (Careful) is always
# correct-high and item3 (Normal) is always correct-low (see codebook note
# in simulate_panel_test.R).
item1_correct_high <- c(TRUE, TRUE, FALSE, FALSE, FALSE)

is_correct <- function(x, correct_high) {
  code <- suppressWarnings(as.numeric(x))
  # %in% resolves an NA lookup to FALSE rather than NA, so it must be
  # restored explicitly -- otherwise missing items are silently scored as
  # "incorrect" instead of propagating as missing.
  result <- if (correct_high) code %in% c(4, 5) else code %in% c(1, 2)
  result[is.na(code)] <- NA
  result
}

scenario_items <- map(1:5, function(s) {
  row <- scenario_cols[s, ]
  tibble(
    !!paste0("scenario", s, "_item1")        := decode(panel_data[[row$item1]], likelihood_labels, ordered = TRUE),
    !!paste0("scenario", s, "_item2_careful") := decode(panel_data[[row$item2]], likelihood_labels, ordered = TRUE),
    !!paste0("scenario", s, "_item3_normal")  := decode(panel_data[[row$item3]], likelihood_labels, ordered = TRUE),
    !!paste0("scenario", s, "_enough_info")   := decode(panel_data[[paste0("Scenario_", s, "_enough_info")]],
                                                          enough_info_labels, ordered = TRUE, level_order = enough_info_order)
  )
}) |> bind_cols()

# primary_correct_count: 0-15 correct across the 5 scenarios x 3 behaviour
# items (item1 + Careful + Normal). Computed on the raw numeric codes.
# do.call(cbind, ...) stacks the 5 (n x 3) matrices side by side into one
# (n x 15) matrix -- using sapply() directly here would instead return a
# 3-D array (n x 3 x 5), which rowSums() would silently mis-sum.
primary_items <- do.call(cbind, lapply(1:5, function(s) {
  row <- scenario_cols[s, ]
  cbind(
    is_correct(panel_data[[row$item1]], item1_correct_high[s]),
    is_correct(panel_data[[row$item2]], TRUE),
    is_correct(panel_data[[row$item3]], FALSE)
  )
}))
# Rows where all 15 primary-outcome items are missing (survey dropout
# before reaching the outcome questions) should stay NA rather than
# collapse to 0 via na.rm; rows with only some items missing (genuine
# partial non-response) still sum over what's available. In practice none
# of the four raw files contain fully-missing rows (see header note), but
# this is kept in case genuine partial non-response exists.
primary_all_missing   <- rowSums(is.na(primary_items)) == ncol(primary_items)
primary_correct_count <- ifelse(primary_all_missing, NA_real_, rowSums(primary_items, na.rm = TRUE))
primary_threshold_12  <- primary_correct_count >= 12

sufficient_info_count <- rowSums(
  sapply(1:5, function(s) {
    code <- suppressWarnings(as.numeric(panel_data[[paste0("Scenario_", s, "_enough_info")]]))
    code %in% c(3, 4)
  }),
  na.rm = TRUE
)

## ---- Other derived outcomes ------------------------------------------------

infant_correct <- with(panel_data,
  Infant.1 == "1" & Infant.2 == "1" & Infant.3 == "0" & Infant.4 == "1" & Infant.5 == "0"
)
high_risk_correct <- with(panel_data,
  High_risk.1 == "1" & High_risk.2 == "1" & High_risk.6 == "1" &
    High_risk.3 == "0" & High_risk.4 == "1" & High_risk.5 == "0"
)
main_aim_correct <- panel_data$Main_aim_advice == "2"

## ---- Assemble analysis-ready tibble ---------------------------------------

panel_data <- panel_data |>
  transmute(
    participant_id = paste0("P", `$submission_id`),
    arm            = arm,
    submitted_at   = `$created`,
    answer_time_ms = `$answer_time_ms`,

    # No raw-data equivalent (see header note) -- everyone treated as
    # having met the inclusion criteria.
    included = Included == "1",

    gender    = decode(Gender, gender_labels),
    age_group = decode(Age, age_labels, ordered = TRUE),
    region    = decode(Municipality, region_labels),
    education = decode(Education, education_labels, ordered = TRUE),
    research_participation_experience = decode(research_participation_experience, research_exp_labels),

    # Employment.1 is the "currently employed" checkbox, coded 0 = not
    # selected / 1 = selected (NOT the 1 = Yes / 2 = No scale used by
    # yes_no_labels elsewhere). It also gates the skip logic below.
    employed                    = factor(ifelse(Employment.1 == "1", "Yes", "No"), levels = c("No", "Yes")),
    contact_high_risk_patients  = decode(Close_contact_patients, yes_no_dontknow_labels),
    can_work_from_home          = decode(Home_office, yes_no_dontknow_labels),
    baseline_behaviour          = decode(Last_time_sick, last_time_sick_labels),
    last_time_home_office       = decode(Last_time_home_office, yes_no_labels),

    reasons_home_office_too_sick = decode(`Reasons_home_office.too_sick`, extent_labels, ordered = TRUE),
    reasons_home_office_infect   = decode(`Reasons_home_office.infect`, extent_labels, ordered = TRUE),
    reasons_home_office_stigma   = decode(`Reasons_home_office.stigma`, extent_labels, ordered = TRUE),
    reasons_home_office_anyway   = decode(`Reasons_home_office.home_office_anyway`, extent_labels, ordered = TRUE),

    sick_note                  = decode(Sick_note, yes_no_labels),
    reasons_sick_note_too_sick = decode(`Reasons_sick_note.too_sick`, extent_labels, ordered = TRUE),
    reasons_sick_note_infect   = decode(`Reasons_sick_note.infect`, extent_labels, ordered = TRUE),
    reasons_sick_note_stigma   = decode(`Reasons_sick_note.stigma`, extent_labels, ordered = TRUE),

    health_literacy = decode(Health_literacy, health_lit_labels, ordered = TRUE),

    main_aim_advice         = decode(Main_aim_advice, main_aim_labels),
    confident_understanding = decode(Confident_understanding, confidence_labels, ordered = TRUE),

    # Infant.*/High_risk.*: binary multi-select checkboxes -> logical.
    infant_1 = Infant.1 == "1", infant_2 = Infant.2 == "1", infant_3 = Infant.3 == "1",
    infant_4 = Infant.4 == "1", infant_5 = Infant.5 == "1",
    high_risk_1 = High_risk.1 == "1", high_risk_2 = High_risk.2 == "1",
    high_risk_3 = High_risk.3 == "1", high_risk_4 = High_risk.4 == "1",
    high_risk_5 = High_risk.5 == "1", high_risk_6 = High_risk.6 == "1",

    trustworthiness      = decode(Trustworthy, trust_labels, ordered = TRUE),
    perceived_usefulness = decode(Seek_info, agreement_labels, ordered = TRUE),
    intention_to_share   = decode(Share_info, agreement_labels, ordered = TRUE),

    # Placeholder: no attention-check item exists in the raw export.
    # Everyone is treated as passing until a real item is added/found.
    attention_check_pass = TRUE
  ) |>
  bind_cols(scenario_items) |>
  bind_cols(tibble(
    # --- Derived outcomes (per codebook "Derived outcomes" sheet) ---
    primary_correct_count = primary_correct_count,
    primary_threshold_12  = primary_threshold_12,
    infant_correct        = infant_correct,
    high_risk_correct     = high_risk_correct,
    main_aim_correct      = main_aim_correct,
    sufficient_info_count = sufficient_info_count
  ))

# Complete-case version: one row per actual survey completion, no
# placeholder rows for the randomised-but-missing participants added below.
# Kept as a separate artifact for reference/debugging -- the rest of the
# pipeline reads "panel_test.rds" (below), not this file.
saveRDS(panel_data, "panel_cc.rds")

## ---- Add back the randomised-but-missing participants ---------------------
# The four raw export files contain only completed submissions (see header
# note) -- respondents who were randomised but never answered (or dropped
# out before) the survey simply aren't rows in them at all. Every downstream
# script, however, expects `panel` to contain ALL randomised participants
# per arm (that's what the CONSORT "Missing data" box and Table 1's
# explicit "Missing" row are built from -- see simulate_panel_test.R's
# convention, replicated here: a missing-outcome row has every field blank
# except participant_id/arm/included). Pad each arm back up to its true
# n_allocated_by_arm with such placeholder rows so per-arm row counts in
# the analysis dataset match the real randomisation counts.
n_completers_by_arm <- table(panel_data$arm)[names(n_allocated_by_arm)]
n_missing_by_arm     <- n_allocated_by_arm - as.integer(n_completers_by_arm)
stopifnot(all(n_missing_by_arm >= 0))   # completers can't exceed those allocated
n_missing_by_arm

missing_rows <- map_dfr(names(n_missing_by_arm), function(a) {
  n_miss <- n_missing_by_arm[[a]]
  if (n_miss == 0) return(tibble())
  tibble(
    # Clearly-synthetic ID: these participants have no raw-data row at all
    # (not even a submission id), unlike real participant_ids ("P<digits>").
    participant_id = paste0("MISSING_", a, "_", seq_len(n_miss)),
    arm            = a,
    # Randomised (hence counted in n_allocated) implies they met the
    # inclusion criteria -- everything else about them is unobserved.
    included       = TRUE
  )
})

# bind_rows() fills every column not present in missing_rows (i.e.
# everything except participant_id/arm/included) with NA of the matching
# type, which is exactly the "missing every field" convention above.
panel_data <- bind_rows(panel_data, missing_rows)
stopifnot(all(as.integer(table(panel_data$arm)[names(n_allocated_by_arm)]) == n_allocated_by_arm))

panel_data |> count(arm)

# Saved under the SAME filename the rest of the pipeline already expects
# (see the "NB" note at the top of descreptive_analysis.R: every downstream
# script reads `panel <- readRDS("panel_test.rds")`) so no other script
# needs to change to run on this real, blinded data. This overwrites the
# simulated test fixture -- re-run simulate_panel_test.R + the previous
# version of this script (or restore from version control) to get it back.
saveRDS(panel_data, "panel_data.rds")

# Real per-arm "allocated" counts (via the blinded file -> arm mapping
# above), for the CONSORT flow chart in descreptive_analysis.R. Now also
# directly verifiable as table(panel_test$arm), but kept as its own small
# file since that's where the "true" number originates (and it's cheap
# insurance against panel_test.rds being rebuilt differently later).
saveRDS(n_allocated_by_arm, "panel_enrolment.rds")
