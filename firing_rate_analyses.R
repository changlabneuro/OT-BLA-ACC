# Analysis of Neural Activity Data for Oxytocin Effects on Decision-Making
# This script analyzes multi-unit activity (MUA) data following oxytocin or saline injections into the BLA

# Load required libraries ------------------------------------------------
library(tidyverse)  # Includes dplyr, ggplot2, tidyr for data manipulation and plotting
library(lme4)       # For linear mixed effects models
library(lmerTest)   # For p-values in mixed models
library(emmeans)    # For post-hoc tests

# Import Data -------------------------------------------------------------
# Define path to data files and high prosocial days
# Replace these paths with your repository structure when sharing
DATA_PATH <- "PATH_TO_YOUR_REPOSITORY/data"

# Load list of days classified as "high prosocial" (by preference index)
high_combo_prefindex_days <- read.csv(file.path(DATA_PATH, "high_prosocial_days.csv"))
high_combo_prefindex_days <- high_combo_prefindex_days$days

# Load preprocessed MUA data
targacq_mua <- read.csv(file.path(DATA_PATH, "firing_rates_mua.csv"), header = TRUE)

# Count channels per day/monkey for reference
day_channel_counts <- targacq_mua %>%
  filter(drugs != 'unspecified') %>%
  distinct(days, monkeys) %>%
  group_by(days, monkeys) %>%
  summarise(channel_n = n(), .groups = 'drop')

# Import normalized MUA data files
FOLDER_PATH <- file.path(DATA_PATH, "mua/normalized_to_cueon/targacq_nas_corrected")
file_paths <- list.files(path = FOLDER_PATH, full.names = TRUE)

days_data <- targacq_mua %>%
  filter(drugs != 'unspecified') %>%
  distinct(days)
days <- days_data$days
mua_data <- NULL

# Load and combine data from all days
for (d in 1:length(days)){
  cur_day <- days[d]
  print(paste("Processing", cur_day))
  
  temp_data <- read.csv(file = sprintf("%s/normtocueon_targacq_data_%s.csv", FOLDER_PATH, cur_day))
  mua_data <- rbind(mua_data, temp_data)
}

# Prepare data for analysis
mua_data <- mua_data %>% 
  mutate(state = ifelse(days %in% high_combo_prefindex_days, 'high', 'low')) %>%
  unite(sites, c(days, channel), remove = FALSE)


# Count total sites per brain region for reference
site_counts <- mua_data %>%
  filter(trialtypes == 'choice') %>%
  group_by(region, monkeys) %>%
  summarise(site_count = n_distinct(sites), .groups = 'drop')

# Import significant sites data -------------------------------------------
# Significant sites were previously identified based on activity changes from baseline period
sig_sites_data <- read.csv(file.path(DATA_PATH, 'significant_mua_sites/targacq_-1to1_sig01.csv'))
significant_sites <- unique(sig_sites_data$sites)

# Get lists of significant sites by region
bla_sig_sites <- sig_sites_data %>%
  filter(region == 'bla') %>%
  pull(sites) %>%
  unique()  

acc_sig_sites <- sig_sites_data %>%
  filter(region == 'acc') %>%
  pull(sites) %>%
  unique()

# Direction of overall firing rate effects --------------------------------
mean_level_effects <- sig_sites_data %>%
  filter(outcomes == 'none', administration == 'post') %>%
  mutate(state = ifelse(days %in% high_combo_prefindex_days, 'high', 'low')) %>%
  filter(case_when(administration == 'pre' ~ abs_trial_num <= 200,
                   administration == 'post' ~ abs_trial_num <= 400)) %>%
  group_by(monkeys, region, outcomes, state, administration, drugs) %>%
  summarise(mean_mua = mean(norm_targacq_psth), .groups = 'drop') %>%
  pivot_wider(names_from = drugs, values_from = mean_mua) %>%
  mutate(sal_min_ot = saline - oxytocin) %>%
  arrange(region, state, monkeys)

# Mixed Effects Model Analysis -------------------------------------------

# Prepare data for mixed effects model 
mean_targacq_mua <- sig_sites_data %>%
  mutate(state = ifelse(days %in% high_combo_prefindex_days, 'high', 'low')) %>%
  filter(case_when(administration=='pre' ~ abs_trial_num <= 200,
                   administration=='post' ~ abs_trial_num <= 400)) %>%
  filter(outcomes == 'other' | outcomes == 'none', administration == 'post', region == 'bla') %>%
  mutate(outcomes = factor(outcomes)) %>%
  mutate(drugs = factor(drugs, levels = c("oxytocin", "saline"))) %>%
  mutate(state = factor(state, levels = c("high", "low")))

# Run full mixed effects model
mixed_model_outcomes <- lmer(norm_targacq_psth ~ drugs * time_number * outcomes * state + 
                               (1|days) + (1|monkeys), 
                             data = mean_targacq_mua)

# Get model summary
summary(mixed_model_outcomes)

# Get estimated marginal means for key combinations
emmeans_comb <- emmeans(mixed_model_outcomes, ~ drugs | state | outcomes)

# Run pairwise comparisons 
pairs(emmeans_comb, adjust = "tukey")


