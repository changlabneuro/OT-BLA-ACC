
# This script analyzes behavioral data from experiments testing how oxytocin infusion in the BLA affects
# decision-making preferences (compared to saline infusion) and performance metrics in monkeys

# Load required libraries ------------------------------------------------

library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(lme4)       # For linear mixed effects models
library(emmeans)    # For estimated marginal means

# Import Data -------------------------------------------------------------

# Load behavioral data, excluding trials where monkeys errored out
all_labels <- read.csv("PATH_TO_YOUR_REPOSITORY/behavioral_data.csv")

# Load behavioral data, all trials including errors 
all_labels_errors <- read.csv("PATH_TO_YOUR_REPOSITORY/behavioral_data_with_errors.csv")

# Load reaction time data
rt_labels <- read.csv("PATH_TO_YOUR_REPOSITORY/reaction_time_data.csv")

# Load state data (assigning states to each day)
high_combo_prefindex_days <- read.csv("PATH_TO_YOUR_REPOSITORY/high_combo_prefindex_days.csv")


# SECTION 1: Binned Preference Index Analysis -----------------------------

# Calculate preference indices across binned trials to examine how preference
# changes throughout the experimental session

# Separately bin pre/baseline and post-infusion trials
bins = seq(0, 200, by = 50)
binned_pre_trials <- all_labels %>%
  filter(administration == 'pre', abs_trial_num <= 250) %>% 
  mutate(bins = cut(abs_trial_num, breaks = bins, include.lowest = TRUE))

bins = seq(0, 400, by = 50)
binned_post_trials <- all_labels %>% 
  filter(administration == 'post', abs_trial_num <= 450) %>%
  filter(trialtypes != 'cued', drugs != 'unspecified') %>%
  mutate(bins = cut(abs_trial_num, breaks = bins, include.lowest = TRUE)) 

binned_data <- rbind(binned_pre_trials, binned_post_trials) 


# SECTION 2: Linear Mixed Effects Model for Preference Index --------------

# Analyze how preference index is affected by drug infusion, baseline state,
# and administration period (pre/post) while accounting for random effects

session_pref_index <- binned_data %>%
  mutate(state = ifelse(days %in% high_combo_prefindex_days, 'high', 'low')) %>%
  filter(outcomes != "errors", drugs != "unspecified", trialtypes != "cued", !is.na(bins)) %>%
  select(administration, contexts, drugs, days, outcomes, trialtypes, abs_trial_num, bins, state, monkeys) %>%
  # Count trials per bin and drug condition
  group_by(bins, drugs, administration, state, monkeys) %>%
  mutate(sample_size = n()) %>%
  # Count outcomes by type
  group_by(drugs, days, contexts, bins, administration, outcomes, state, monkeys) %>%
  mutate(outcome_count = n()) %>%
  mutate(self_count = ifelse(outcomes == 'self', outcome_count, NA)) %>%
  mutate(both_count = ifelse(outcomes == 'both', outcome_count, NA)) %>%
  mutate(other_count = ifelse(outcomes == 'other', outcome_count, NA)) %>%
  mutate(none_count = ifelse(outcomes == 'none', outcome_count, NA)) %>%
  distinct(drugs, days, contexts, administration, bins, self_count, both_count, other_count, none_count, sample_size, state, monkeys) %>%
  # Summarize counts by bin
  group_by(drugs, days, contexts, administration, bins, state, monkeys) %>%
  reframe(sample_size, bins, days, contexts, 
          self_count = sum(self_count, na.rm = TRUE),
          both_count = sum(both_count, na.rm = TRUE), 
          other_count = sum(other_count, na.rm = TRUE),
          none_count = sum(none_count, na.rm = TRUE)) %>%
  distinct(drugs, days, administration, self_count, both_count, other_count, none_count, bins, sample_size, state, monkeys) %>%
  # Calculate preference indices
  mutate(sb_index = (both_count - self_count)/(both_count + self_count)) %>%
  mutate(on_index = (other_count - none_count)/(other_count + none_count)) %>%
  mutate(preference_index = ifelse(is.na(sb_index), on_index, sb_index), 
         contexts = ifelse(is.na(sb_index), 'othernone', 'selfboth')) %>%
  # Factor variables for model
  mutate(drugs = factor(drugs, levels = c("oxytocin", "saline"))) %>%
  mutate(state = factor(state, levels = c("high", "low"))) %>%
  mutate(administration = factor(administration, levels = c("post", "pre"))) %>%
  filter(contexts == 'selfboth')  # Focus on self-both context trials

# Run linear mixed effects model
mixed_model_outcomes <- lmer(preference_index ~ drugs * state * administration * contexts + 
                               (1|days) + (1|monkeys), 
                             data = session_pref_index)

# Summary of the mixed model
summary(mixed_model_outcomes)

# Get estimated marginal means for combinations of drug, state, and outcome
emmeans_comb <- emmeans(mixed_model_outcomes, ~ drugs | state | administration | contexts)

# Pairwise comparisons with Tukey adjustment
pairs(emmeans_comb, adjust = "tukey")



# SECTION 3: Binned Preference Index Statistical Comparisons --------------

# Compare preference indices between drug conditions within each bin using
# non-parametric Wilcoxon rank sum tests

# Prepare data
session_pref_index <- binned_data %>%
  mutate(state = ifelse(days %in% high_combo_prefindex_days, 'high', 'low')) %>%
  filter(outcomes != "errors", drugs != "unspecified", trialtypes != "cued", !is.na(bins)) %>%
  select(administration, contexts, drugs, days, outcomes, trialtypes, abs_trial_num, bins, state) %>%
  # Similar data processing as above but focused on bin comparisons
  group_by(bins, drugs, administration, state) %>%
  mutate(sample_size = n()) %>%
  group_by(drugs, days, contexts, bins, administration, outcomes, state) %>%
  mutate(outcome_count = n()) %>%
  mutate(self_count = ifelse(outcomes == 'self', outcome_count, NA)) %>%
  mutate(both_count = ifelse(outcomes == 'both', outcome_count, NA)) %>%
  mutate(other_count = ifelse(outcomes == 'other', outcome_count, NA)) %>%
  mutate(none_count = ifelse(outcomes == 'none', outcome_count, NA)) %>%
  distinct(drugs, days, contexts, administration, bins, self_count, both_count, other_count, none_count, sample_size, state) %>%
  group_by(drugs, days, contexts, administration, bins, state) %>%
  reframe(sample_size, bins, days, contexts, 
          self_count = sum(self_count, na.rm = TRUE),
          both_count = sum(both_count, na.rm = TRUE), 
          other_count = sum(other_count, na.rm = TRUE),
          none_count = sum(none_count, na.rm = TRUE)) %>%
  distinct(drugs, days, administration, self_count, both_count, other_count, none_count, bins, sample_size, state) %>%
  mutate(sb_index = (both_count - self_count)/(both_count + self_count)) %>%
  mutate(on_index = (other_count - none_count)/(other_count + none_count)) %>%
  mutate(preference_index = ifelse(is.na(sb_index), on_index, sb_index), 
         contexts = ifelse(is.na(sb_index), 'othernone', 'selfboth')) %>%
  filter(contexts == 'othernone')  # Focus on other-none context trials


# SECTION 4: Overall Preference Index Summary -----------------------------

# Calculate overall preference indices summarized by drug, state, and monkey

pref_index_summary <- all_labels %>% 
  filter(outcomes != "errors", drugs != "unspecified", trialtypes != "cued") %>%
  filter((administration == 'pre' & abs_trial_num <= 200) | 
           (administration == 'post' & abs_trial_num <= 400)) %>%
  mutate(state = ifelse(days %in% high_combo_prefindex_days, 'high', 'low')) %>%
  ungroup() %>%
  select(administration, contexts, drugs, days, outcomes, trialtypes, trial_number, state, monkeys) %>%
  group_by(drugs, administration, contexts, state, monkeys) %>%
  mutate(sample_size = n()/49) %>%  # Adjust divisor based on your experimental design
  group_by(drugs, days, contexts, administration, outcomes, state) %>%
  mutate(outcome_count = n()) %>%
  mutate(self_count = ifelse(outcomes == 'self', outcome_count, NA)) %>%
  mutate(both_count = ifelse(outcomes == 'both', outcome_count, NA)) %>%
  mutate(other_count = ifelse(outcomes == 'other', outcome_count, NA)) %>%
  mutate(none_count = ifelse(outcomes == 'none', outcome_count, NA)) %>%
  distinct(drugs, days, contexts, administration, self_count, both_count, other_count, none_count, sample_size, state, monkeys) %>%
  group_by(drugs, days, contexts, administration, state, monkeys) %>%
  summarise(sample_size = first(sample_size), 
            days = first(days), 
            contexts = first(contexts),
            self_count = sum(self_count, na.rm = TRUE),
            both_count = sum(both_count, na.rm = TRUE), 
            other_count = sum(other_count, na.rm = TRUE),
            none_count = sum(none_count, na.rm = TRUE),
            .groups = 'drop') %>%
  distinct(drugs, days, administration, self_count, both_count, other_count, none_count, sample_size, monkeys) %>%
  mutate(sb_index = (both_count - self_count)/(both_count + self_count)) %>%
  mutate(on_index = (other_count - none_count)/(other_count + none_count)) %>%
  group_by(drugs, administration, state, days, monkeys) %>%
  mutate(preference_index = ifelse(is.na(sb_index), on_index, sb_index)) %>%
  group_by(administration, drugs, contexts, state) %>%
  mutate(sd_index = sd(preference_index, na.rm = TRUE)) %>%
  mutate(sem_index = sd_index/sqrt(sample_size))

# Statistical test example - comparing pre vs post for a specific condition
# Input parameters for specific comparison
drug <- 'oxytocin'
cur_state <- 'low'
context <- 'selfboth'

pre_data <- pref_index_summary %>% 
  filter(drugs == drug, state == cur_state, contexts == context, administration == 'pre')
post_data <- pref_index_summary %>% 
  filter(drugs == drug, state == cur_state, contexts == context, administration == 'post')

wilcoxon_test <- wilcox.test(pre_data$preference_index, post_data$preference_index, paired = FALSE)
print(wilcoxon_test)


# SECTION 5: Error Rate Analysis Across Session ---------------------------

# Analyze how error rates change across the session by bin, drug, and state

# Prepare binned data including error trials
bins = seq(0, 200, by = 50)
binned_pre_trials <- all_labels_errors %>%
  filter(administration == 'pre' & abs_trial_num <= 200) %>% 
  mutate(bins = cut(abs_trial_num, breaks = bins, include.lowest = TRUE)) 

bins = seq(0, 400, by = 50)
binned_post_trials <- all_labels_errors %>% 
  filter(administration == 'post' & abs_trial_num <= 400) %>%
  mutate(bins = cut(abs_trial_num, breaks = bins, include.lowest = TRUE)) 

binned_data <- rbind(binned_pre_trials, binned_post_trials) 

# Calculate error rates by bin
errors_across_trials <- binned_data %>%
  filter(trialtypes == 'choice', drugs != 'unspecified') %>%
  mutate(state = ifelse(days %in% high_combo_prefindex_days, 'high', 'low')) %>%
  group_by(administration, contexts, trialtypes, drugs, state, bins, days, monkeys) %>%
  mutate(trial_count = n()) %>%
  mutate(error_count = sum(outcomes == 'errors')) %>% 
  distinct(administration, contexts, trialtypes, drugs, state, bins, days, error_count, trial_count, monkeys) %>%
  mutate(day_error_prop = error_count/trial_count) %>%
  group_by(administration, state, drugs, trialtypes, contexts, bins) %>%
  summarise(mean_error_prop = mean(day_error_prop, na.rm = TRUE), 
            n_n = n(), 
            sd = sd(day_error_prop, na.rm = TRUE), 
            sem = sd/sqrt(n_n),
            .groups = 'drop')

# Mixed effects model for error rates
lme_data <- errors_across_trials %>% 
  filter(contexts == 'othernone') %>%
  mutate(drugs = factor(drugs, levels = c("oxytocin", "saline"))) %>%
  mutate(state = factor(state, levels = c("high", "low"))) %>%
  mutate(administration = factor(administration, levels = c("post", "pre"))) 

mixed_model_errors <- lmer(day_error_prop ~ drugs * state * administration + (1|days) + (1|monkeys), 
                           data = lme_data)

# Summary of the mixed model
summary(mixed_model_errors)

# Get estimated marginal means
emmeans_errors <- emmeans(mixed_model_errors, ~ drugs | state | administration)
pairs(emmeans_errors, adjust = "tukey")



# SECTION 6: Reaction Time Analysis ---------------------------------------

# Analyze how reaction times are affected by drug, state, and outcome type

# Process the reaction time data
rt_label_no_errors <- rt_labels %>%
  filter(outcomes != 'errors') %>%
  mutate(state = ifelse(days %in% high_combo_prefindex_days, 'high', 'low')) %>%
  separate_rows(trials, sep = "_") %>% 
  filter(trials != 'trial', trials != '') %>%
  mutate(trial_number = as.numeric(trials))

rt_trial_data <- rt_label_no_errors %>%
  mutate(state = ifelse(days %in% high_combo_prefindex_days, 'high', 'low')) %>%
  distinct(administration, outcomes, drugs, days, blocks, trial_number, monkeys, contexts, trialtypes, state, reaction_time) %>%
  filter(drugs !='unspecified') %>%
  # Clean up context labels
  mutate(contexts = ifelse(outcomes == 'other' | outcomes == 'none', 'othernone', 
                           ifelse(outcomes == 'self' | outcomes == 'both', 'selfboth', contexts))) %>%
  # Standardize context naming
  mutate(contexts = ifelse(contexts %in% c('context__both', 'context__self'), 'selfboth', contexts)) %>% 
  mutate(contexts = ifelse(contexts %in% c('context__other', 'context__none'), 'othernone', contexts)) %>%
  mutate(abs_trial_num = trial_number)

# Calculate absolute trial number
for (i in 2:nrow(rt_trial_data)) {
  if (rt_trial_data[i,'administration'] == rt_trial_data[i-1,'administration'] & 
      (rt_trial_data[i,'blocks'] != rt_trial_data[i-1,'blocks'])) {
    trial_num_point = rt_trial_data[i-1,'abs_trial_num']
    rt_trial_data[i,'abs_trial_num'] = rt_trial_data[i,'trial_number'] + trial_num_point + 1
    
  } else if (rt_trial_data[i,'administration'] != rt_trial_data[i-1,'administration'] & 
             (rt_trial_data[i,'blocks'] != rt_trial_data[i-1,'blocks'])) {
    rt_trial_data[i,'abs_trial_num'] = rt_trial_data[i,'trial_number']
    
  } else if (rt_trial_data[i,'administration'] == rt_trial_data[i-1,'administration'] & 
             (rt_trial_data[i,'blocks'] == rt_trial_data[i-1,'blocks'])) {
    rt_trial_data[i,'abs_trial_num'] = (rt_trial_data[i-1,'abs_trial_num'] + 
                                          (rt_trial_data[i,'trial_number'] - rt_trial_data[i-1,'trial_number']))
  }
}

# Bin reaction time data
bins = seq(0, 200, by = 50)
binned_pre_trials <- rt_trial_data %>%
  filter(administration == 'pre', abs_trial_num <= 200) %>% 
  mutate(bins = cut(abs_trial_num, breaks = bins, include.lowest = TRUE)) 

bins = seq(0, 400, by = 50)
binned_post_trials <- rt_trial_data %>% 
  filter(administration == 'post', abs_trial_num <= 400) %>%
  mutate(bins = cut(abs_trial_num, breaks = bins, include.lowest = TRUE)) %>%
  group_by(days, bins, contexts, administration, drugs, trialtypes) %>%
  mutate(bin_trial_count = n()) %>%
  select(-bin_trial_count)

# Combine pre and post trial data
rt_binned_data <- rbind(binned_pre_trials, binned_post_trials) 

# Mixed effects model for reaction time
rt_lme_data <- rt_binned_data %>% 
  filter(outcomes == 'other' | outcomes == 'none') %>%
  mutate(drugs = factor(drugs, levels = c("oxytocin", "saline"))) %>%
  mutate(state = factor(state, levels = c("high", "low"))) %>%
  mutate(administration = factor(administration, levels = c("post", "pre")))

# Run the linear mixed effects model for reaction time
mixed_model_rt <- lmer(reaction_time ~ drugs * state * administration * outcomes + 
                         (1|days) + (1|monkeys), 
                       data = rt_lme_data)

# Summary of the model
summary(mixed_model_rt)

# Get estimated marginal means
emmeans_rt <- emmeans(mixed_model_rt, ~ drugs | state | administration | outcomes)
pairs(emmeans_rt, adjust = "tukey")
