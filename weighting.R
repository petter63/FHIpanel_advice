# Post-stratification (raking) weights, computed SEPARATELY within each
# randomised arm, targeting national population margins for
# age x gender (joint) and education level (marginal).
#
# Population sources (SSB):
#  - "Alder_07459_20260921-122645.csv": Befolkning etter region, alder,
#    kjonn og aar (2026), whole country ("0 Hele landet").
#  - "Utdanningsnivaa for personer 16 aar og eldre.csv": Utdanningsniva,
#    2020 and 2025, whole country, ages 16+.
#
# Cleaning decisions confirmed with the user (2026-09-21):
#  - Education: the population table only has 4 levels (no split of
#    university/college into short/long cycle), so the survey's 5-level
#    `education` factor is collapsed to 4 levels to match, by merging
#    "University/college, 4 years or less" and "more than 4 years" into a
#    single "University/college" category for RAKING PURPOSES ONLY. The
#    original 5-level `education` is kept as-is for descriptive reporting.
#  - Gender: the population table has no "Other" / "Prefer not to say"
#    category. Respondents in these groups get weight = 1 (no
#    post-stratification adjustment), since there is no population
#    benchmark to rake them against.
#  - Age: SSB's "15-19" 5-year bin is mapped onto the survey's "16-24"
#    group (it also contains 15-year-olds, who are outside the survey's
#    sampling frame) -- a minor approximation.
#  - Respondents with a missing age_group, gender, or education value are
#    excluded from raking (weight = NA); flagged for review before
#    analysis, since substantial missingness on these variables would bias
#    the weighted estimates.

library(tidyverse)
library(survey)

pop_dir <- "C:/Users/peel/OneDrive - Folkehelseinstituttet/Studier/Panel_smittevernraad"

panel_test <- readRDS("panel_test.rds")

## ---- 1. Build population margins ------------------------------------------

age_raw <- read_delim(
  file.path(pop_dir, "Alder_07459_20260921-122645.csv"),
  delim = ";", skip = 1, locale = locale(encoding = "Latin1"),
  col_types = cols(.default = "c")
)
names(age_raw) <- c("region", "age_raw", "gender_raw", "n")
age_raw <- age_raw |> mutate(n = as.numeric(n))

age_bin_to_group <- function(x) {
  case_when(
    x %in% c("15-19 \u00e5r", "20-24 \u00e5r") ~ "16-24",
    x %in% c("25-29 \u00e5r", "30-34 \u00e5r") ~ "25-34",
    x %in% c("35-39 \u00e5r", "40-44 \u00e5r") ~ "35-44",
    x %in% c("45-49 \u00e5r", "50-54 \u00e5r") ~ "45-54",
    x %in% c("55-59 \u00e5r", "60-64 \u00e5r") ~ "55-64",
    x %in% c("65-69 \u00e5r", "70-74 \u00e5r", "75-79 \u00e5r", "80-84 \u00e5r",
             "85-89 \u00e5r", "90-94 \u00e5r", "95-99 \u00e5r",
             "100 \u00e5r eller eldre") ~ "65 or older",
    TRUE ~ NA_character_  # below the survey's sampling frame (0-14)
  )
}

age_group_levels <- c("16-24", "25-34", "35-44", "45-54", "55-64", "65 or older")

age_gender_pop <- age_raw |>
  filter(region == "0 Hele landet") |>
  mutate(
    age_group = factor(age_bin_to_group(age_raw), levels = age_group_levels, ordered = TRUE),
    gender    = recode(gender_raw, "Kvinner" = "Female", "Menn" = "Male")
  ) |>
  filter(!is.na(age_group))

# Marginal (not joint) age and gender population shares. With ~250
# respondents per arm split across 6 age groups x 2 genders, some
# age x gender cells are empty within an arm, which breaks joint raking
# (rake()/postStratify() cannot rake a population cell that has zero
# sample support). Raking on the two separate one-way margins avoids
# this and is standard practice at this sample size; it matches the
# age and gender marginal distributions but not necessarily their
# cross-tabulation within each arm.
age_pop <- age_gender_pop |>
  summarise(n = sum(n), .by = age_group) |>
  mutate(p = n / sum(n)) |>
  select(age_group, p)

gender_pop <- age_gender_pop |>
  summarise(n = sum(n), .by = gender) |>
  mutate(p = n / sum(n)) |>
  select(gender, p)

edu_raw <- read_delim(
  file.path(pop_dir, "Utdanningsniv\u00e5 for personer 16 \u00e5r og eldre.csv"),
  delim = ";", locale = locale(encoding = "UTF-8"),
  col_types = cols(.default = "c")
)

education4_levels <- c("Primary school or lower", "Upper secondary school",
                        "Vocational college", "University/college")

edu_pop <- edu_raw |>
  slice(2:6) |>
  transmute(
    category = `...1`,
    n_2025   = as.numeric(gsub("\\s", "", `2025`))
  ) |>
  filter(category != "Utdanningsniv\u00e5, i alt") |>
  mutate(education4 = factor(recode(category,
      "Grunnskoleniv\u00e5"                    = "Primary school or lower",
      "Videreg\u00e5ende skoleniv\u00e5"       = "Upper secondary school",
      "Fagskoleniv\u00e5"                      = "Vocational college",
      "Universitets- og h\u00f8gskoleniv\u00e5" = "University/college"
    ), levels = education4_levels)) |>
  summarise(n = sum(n_2025), .by = education4) |>
  mutate(p = n / sum(n)) |>
  select(education4, p)

## ---- 2. Prepare survey data for raking -------------------------------------

panel_weighting <- panel_test |>
  mutate(
    # Collapse to 4 levels to match the population education margin;
    # `education` (5 levels) is left untouched for descriptive use.
    education4 = fct_collapse(
      education,
      "University/college" = c("University/college, 4 years or less",
                                "University/college, more than 4 years")
    ) |> fct_relevel(education4_levels),
    # No population benchmark for these groups -> excluded from raking,
    # handled separately below (weight = 1).
    no_pop_benchmark = gender %in% c("Other", "Prefer not to say") | is.na(gender),
    raking_complete  = !no_pop_benchmark & !is.na(age_group) & !is.na(education4)
  )

n_excluded <- sum(!panel_weighting$raking_complete)
if (n_excluded > 0) {
  message(n_excluded, " respondent(s) excluded from raking (missing age/gender/",
          "education, or gender = Other/Prefer not to say).")
}

## ---- 3. Rake separately within each arm ------------------------------------

rake_arm <- function(df) {
  # Only drop unused levels on the raking variables themselves -- calling
  # droplevels() on the whole data frame would also strip levels from
  # unrelated ordered factors, making them incompatible to row-bind back
  # together across arms.
  df_complete <- df |>
    filter(raking_complete) |>
    mutate(age_group = droplevels(age_group), gender = droplevels(gender))
  n_arm <- nrow(df_complete)

  # Population margins rescaled to this arm's (complete-case) sample size,
  # so all margin tables sum to the same total, as required by rake().
  age_margin <- age_pop |>
    filter(age_group %in% levels(df_complete$age_group)) |>
    mutate(Freq = p / sum(p) * n_arm) |>
    select(age_group, Freq)

  gender_margin <- gender_pop |>
    filter(gender %in% levels(df_complete$gender)) |>
    mutate(Freq = p / sum(p) * n_arm) |>
    select(gender, Freq)

  education_margin <- edu_pop |>
    mutate(Freq = p * n_arm) |>
    select(education4, Freq)

  design <- svydesign(ids = ~1, data = df_complete, weights = ~1)
  raked  <- rake(
    design,
    sample.margins     = list(~age_group, ~gender, ~education4),
    population.margins = list(age_margin, gender_margin, education_margin)
  )

  df_complete$weight <- weights(raked)
  df_complete
}

panel_weighting <- panel_weighting |>
  group_split(arm) |>
  map(rake_arm) |>
  list_rbind() |>
  # Restore the excluded rows with weight = 1 (no adjustment; see header note)
  bind_rows(panel_weighting |> filter(!raking_complete) |> mutate(weight = 1)) |>
  select(-no_pop_benchmark, -raking_complete)

## ---- 4. Sanity checks -------------------------------------------------------

# Weights should average ~1 within each arm and not be extreme.
panel_weighting |>
  summarise(
    n = n(), mean_weight = mean(weight), sd_weight = sd(weight),
    min_weight = min(weight), max_weight = max(weight),
    .by = arm
  )

saveRDS(panel_weighting, "panel_weighting.rds")
