## Clean the raw survey-export test data (panel_test_raw: a single data
## frame with one row per respondent, an "arm" column, and otherwise the
## same raw variables as the real "Smittevern *.xlsx" files, as produced by
## simulate_panel_test.R) into the analysis-ready tibble structure that
## descreptive_analysis.R expects.
##
## Value labels and derived outcomes below come from
## "data/Smittevern_RCT_codebook and data dictionary.xlsx" (sheets
## "Detailed variable key" and "Derived outcomes"). A few things to flag:
##
##  - CRITICAL (per the codebook's own "Review before analysis" sheet):
##    the real, verified randomised allocation is *not* in the raw export.
##    `arm` here is still the placeholder file->arm mapping from
##    simulate_panel_test.R, not a verified allocation. Replace once the
##    real randomisation log is merged in.
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
##  - attention_check_pass has no raw-data equivalent in the export, so it
##    is set to TRUE for everyone (no exclusions) until a real
##    attention-check item is identified in the survey.

library(tidyverse)

panel_test_raw <- readRDS("panel_test_raw.rds")

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

# Column names/casing are fixed per scenario (see simulate_panel_test.R);
# all three items per scenario share the same "likelihood" response scale.
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
  if (correct_high) code %in% c(4, 5) else code %in% c(1, 2)
}

scenario_items <- map(1:5, function(s) {
  row <- scenario_cols[s, ]
  tibble(
    !!paste0("scenario", s, "_item1")        := decode(panel_test_raw[[row$item1]], likelihood_labels, ordered = TRUE),
    !!paste0("scenario", s, "_item2_careful") := decode(panel_test_raw[[row$item2]], likelihood_labels, ordered = TRUE),
    !!paste0("scenario", s, "_item3_normal")  := decode(panel_test_raw[[row$item3]], likelihood_labels, ordered = TRUE),
    !!paste0("scenario", s, "_enough_info")   := decode(panel_test_raw[[paste0("Scenario_", s, "_enough_info")]],
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
    is_correct(panel_test_raw[[row$item1]], item1_correct_high[s]),
    is_correct(panel_test_raw[[row$item2]], TRUE),
    is_correct(panel_test_raw[[row$item3]], FALSE)
  )
}))
primary_correct_count <- rowSums(primary_items, na.rm = TRUE)
primary_threshold_12  <- primary_correct_count >= 12

sufficient_info_count <- rowSums(
  sapply(1:5, function(s) {
    code <- suppressWarnings(as.numeric(panel_test_raw[[paste0("Scenario_", s, "_enough_info")]]))
    code %in% c(3, 4)
  }),
  na.rm = TRUE
)

## ---- Other derived outcomes ------------------------------------------------

infant_correct <- with(panel_test_raw,
  Infant.1 == "1" & Infant.2 == "1" & Infant.3 == "0" & Infant.4 == "1" & Infant.5 == "0"
)
high_risk_correct <- with(panel_test_raw,
  High_risk.1 == "1" & High_risk.2 == "1" & High_risk.6 == "1" &
    High_risk.3 == "0" & High_risk.4 == "1" & High_risk.5 == "0"
)
main_aim_correct <- panel_test_raw$Main_aim_advice == "2"

## ---- Assemble analysis-ready tibble ---------------------------------------

panel_test <- panel_test_raw |>
  transmute(
    participant_id = paste0("P", `$submission_id`),
    arm            = arm,
    submitted_at   = `$created`,
    answer_time_ms = `$answer_time_ms`,

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

saveRDS(panel_test, "panel_test.rds")
