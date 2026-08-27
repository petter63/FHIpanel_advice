## Simulate test data matching the Panel_smittevernråd trial protocol
## (baseline characteristics, primary outcome, and secondary outcomes)

library(tidyverse)

set.seed(4831)

n <- 200
# Arm coding per trial protocol:
#   V1 (control)            - current formulation
#   V2                      - + sentence on when it's okay to work/participate in activities
#   V3                      - + definitions of key terms
#   V4                      - + sentence AND definitions
arms <- c("V1_control", "V2_sentence", "V3_definitions", "V4_sentence_definitions")

likert5_extent <- c("Not at all", "To a small extent", "To some extent",
                     "To a large extent", "To a very large extent")
baseline_behav_opts <- c("Not relevant", "Went to work", "Worked from home",
                          "Got a doctor's note", "Used a sick day",
                          "Wasn't working at the time")
wfh_reason_cols <- c("wfh_reason_felt_too_sick", "wfh_reason_worried_infecting_others",
                      "wfh_reason_worried_colleague_reaction", "wfh_reason_planned_wfh_anyway")
sickday_reason_cols <- c("sickday_reason_felt_too_sick", "sickday_reason_worried_infecting_others",
                          "sickday_reason_worried_colleague_reaction", "sickday_reason_planned_wfh_anyway")

# Per-scenario answer key: is the "high likelihood" response the correct one for item_a?
scenario_info <- tibble(
  scenario = 1:5,
  item_a_correct_high = c(TRUE, TRUE, FALSE, FALSE, FALSE)
)

# Assumed probability of answering all 3 items correctly per scenario, by arm
p_correct_by_arm <- c(V1_control = 0.55, V2_sentence = 0.68,
                       V3_definitions = 0.65, V4_sentence_definitions = 0.72)

sim_likert_directional <- function(n, correct_high, p_correct) {
  is_correct <- runif(n) < p_correct
  vapply(seq_len(n), function(i) {
    if (is_correct[i]) {
      if (correct_high) sample(4:5, 1) else sample(1:2, 1)
    } else {
      sample(1:5, 1)
    }
  }, numeric(1))
}

# NOTE: this script assumes `df` has already been assembled with baseline
# characteristics, arm assignment, and the 5-scenario primary outcome items
# (s1_item_a ... s5_enough_info) using the helpers above, before the
# secondary-outcome block below is run. See conversation history / analysis
# notebook for the full baseline + primary-outcome simulation step.

# ---- Secondary outcomes ----

has_explanations <- df$arm %in% c("V3_definitions", "V4_sentence_definitions")

# Understanding of the advice: main goal (Table 4)
goal_opts <- c(
  "Protect yourself",
  "Protect others",
  "Protect others that have a higher risk of becoming seriously ill",
  "Prevent spread of infection in the community"
)
p_goal_correct <- ifelse(has_explanations, 0.62, 0.48)
main_goal_advice <- vapply(seq_len(n), function(i) {
  if (runif(1) < p_goal_correct[i]) {
    goal_opts[3]
  } else {
    sample(goal_opts[-3], 1)
  }
}, character(1))

# Confidence in understanding (1 = Very unsure -> 4 = Very sure)
confidence_understanding <- pmin(4, pmax(1, round(rnorm(n, mean = 2.6, sd = 0.9))))

# Understanding of term "infant": multi-select, correct = options 1, 2, 4
infant_true_opts  <- c(1, 2, 4)
infant_false_opts <- c(3, 5)
p_true_hit  <- ifelse(has_explanations, 0.80, 0.55)
p_false_hit <- ifelse(has_explanations, 0.12, 0.30)

infant_selection <- matrix(FALSE, nrow = n, ncol = 5,
                           dimnames = list(NULL, paste0("infant_opt", 1:5)))
for (i in seq_len(n)) {
  infant_selection[i, infant_true_opts]  <- runif(length(infant_true_opts))  < p_true_hit[i]
  infant_selection[i, infant_false_opts] <- runif(length(infant_false_opts)) < p_false_hit[i]
}
infant_selection <- as_tibble(infant_selection)
infant_term_correct <- apply(infant_selection, 1, function(row) {
  setequal(which(row), infant_true_opts)
})

# Understanding of term "people at higher risk": multi-select, correct = 1, 2, 4, 5
risk_true_opts  <- c(1, 2, 4, 5)
risk_false_opts <- c(3, 6)
p_true_hit_r  <- ifelse(has_explanations, 0.78, 0.50)
p_false_hit_r <- ifelse(has_explanations, 0.15, 0.35)

risk_selection <- matrix(FALSE, nrow = n, ncol = 6,
                          dimnames = list(NULL, paste0("risk_group_opt", 1:6)))
for (i in seq_len(n)) {
  risk_selection[i, risk_true_opts]  <- runif(length(risk_true_opts))  < p_true_hit_r[i]
  risk_selection[i, risk_false_opts] <- runif(length(risk_false_opts)) < p_false_hit_r[i]
}
risk_selection <- as_tibble(risk_selection)
risk_group_term_correct <- apply(risk_selection, 1, function(row) {
  setequal(which(row), risk_true_opts)
})

# Trustworthiness (1-10), Perceived usefulness (1-5), Intention to share (1-5)
trustworthiness      <- pmin(10, pmax(1, round(rnorm(n, mean = 7, sd = 1.8))))
perceived_usefulness <- pmin(5,  pmax(1, round(rnorm(n, mean = 3.6, sd = 1.0))))
intention_to_share   <- pmin(5,  pmax(1, round(rnorm(n, mean = 3.4, sd = 1.1))))

df <- df |>
  bind_cols(
    tibble(
      main_goal_advice        = factor(main_goal_advice, levels = goal_opts),
      main_goal_advice_correct = main_goal_advice == goal_opts[3],
      confidence_understanding = ordered(confidence_understanding, levels = 1:4,
                                          labels = c("Very unsure", "Somewhat unsure",
                                                     "Somewhat sure", "Very sure"))
    ),
    infant_selection,
    tibble(infant_term_correct = infant_term_correct),
    risk_selection,
    tibble(
      risk_group_term_correct = risk_group_term_correct,
      trustworthiness         = trustworthiness,
      perceived_usefulness    = perceived_usefulness,
      intention_to_share      = intention_to_share
    )
  )

panel_test <- df
saveRDS(panel_test, "panel_test.rds")
