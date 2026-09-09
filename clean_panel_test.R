## Clean the raw survey-export test data (panel_test_raw: a single data
## frame with one row per respondent, an "arm" column, and otherwise the
## same raw variables as the real "Smittevern *.xlsx" files, as produced by
## simulate_panel_test.R) into the analysis-ready tibble structure that
## descreptive_analysis.R expects.
##
## Decisions made with the user (2026-09-02), all provisional pending real
## documentation:
##  - No codebook exists yet for the raw numeric codes, so categorical
##    variables are NOT relabeled here -- they keep their raw codes as
##    factors (e.g. gender level "1"/"2"/"3"/"4"). Replace
##    `fct_recode()`/`levels<-()` once real value labels are available.
##  - attention_check_pass has no raw-data equivalent in the export, so it
##    is set to TRUE for everyone (no exclusions) until a real
##    attention-check item is identified in the survey.

library(tidyverse)

panel_test_raw <- readRDS("panel_test_raw.rds")

# Comprehension scenario items: the specific item-1 wording (Delay/Home)
# differs by scenario and hasn't been decoded, so columns are pulled in
# their original left-to-right order and renamed generically
# (item1/item2/item3) rather than guessing semantics.
extract_scenario_items <- function(df, s) {
  item_cols <- grep(paste0("^Scenario_", s, "\\."), names(df), value = TRUE)
  stopifnot(length(item_cols) == 3)
  out <- df |>
    select(all_of(item_cols)) |>
    mutate(across(everything(), factor))
  names(out) <- paste0("scenario", s, "_item", seq_along(item_cols))
  info_col <- paste0("Scenario_", s, "_enough_info")
  out[[paste0("scenario", s, "_enough_info")]] <- factor(df[[info_col]])
  out
}

scenario_items <- map(1:5, ~ extract_scenario_items(panel_test_raw, .x)) |> bind_cols()

panel_test <- panel_test_raw |>
  transmute(
    participant_id = paste0("P", `$submission_id`),
    arm            = arm,
    submitted_at   = `$created`,
    answer_time_ms = `$answer_time_ms`,

    # --- Demographics: raw numeric codes kept as coded factors (see
    #     header note -- no codebook yet to translate to text labels) ---
    gender    = factor(Gender),
    age_group = factor(Age),
    region    = factor(Municipality),
    education = factor(Education),
    research_participation_experience = factor(research_participation_experience),

    # --- Employment & work-related skip logic ---
    # Employment.1 is a binary "in paid work" checkbox (raw code
    # "1"/"0"), kept as a factor like the other coded variables so it
    # stacks cleanly alongside them in downstream summary tables. Raw
    # NAs from the survey's skip logic (only asked when employed) are
    # preserved as-is by simply wrapping the raw values in factor().
    employed                    = factor(Employment.1),
    contact_high_risk_patients  = factor(Close_contact_patients),
    can_work_from_home          = factor(Home_office),
    baseline_behaviour          = factor(Last_time_sick),
    last_time_home_office       = factor(Last_time_home_office),

    reasons_home_office_too_sick = factor(`Reasons_home_office.too_sick`),
    reasons_home_office_infect   = factor(`Reasons_home_office.infect`),
    reasons_home_office_stigma   = factor(`Reasons_home_office.stigma`),
    reasons_home_office_anyway   = factor(`Reasons_home_office.home_office_anyway`),

    sick_note                  = factor(Sick_note),
    reasons_sick_note_too_sick = factor(`Reasons_sick_note.too_sick`),
    reasons_sick_note_infect   = factor(`Reasons_sick_note.infect`),
    reasons_sick_note_stigma   = factor(`Reasons_sick_note.stigma`),

    health_literacy = factor(Health_literacy),

    # --- Secondary outcomes ---
    main_aim_advice          = factor(Main_aim_advice),
    confident_understanding  = factor(Confident_understanding),
    infant_1 = factor(Infant.1), infant_2 = factor(Infant.2),
    infant_3 = factor(Infant.3), infant_4 = factor(Infant.4),
    infant_5 = factor(Infant.5),
    high_risk_1 = factor(High_risk.1), high_risk_2 = factor(High_risk.2),
    high_risk_3 = factor(High_risk.3), high_risk_4 = factor(High_risk.4),
    high_risk_5 = factor(High_risk.5), high_risk_6 = factor(High_risk.6),
    trustworthiness      = factor(Trustworthy),
    perceived_usefulness = factor(Seek_info),
    intention_to_share   = factor(Share_info),

    # Placeholder: no attention-check item exists in the raw export.
    # Everyone is treated as passing until a real item is added/found.
    attention_check_pass = TRUE
  ) |>
  bind_cols(scenario_items)

saveRDS(panel_test, "panel_test.rds")
