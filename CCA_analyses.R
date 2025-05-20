
# This script performs canonical correlation analysis (CCA) between BLA and ACC regions
# to examine how oxytocin affects directional information flow during decision-making

# Load required libraries ------------------------------------------------
library(tidyverse)  # For data manipulation and visualization
library(lme4)       # For linear mixed effects models
library(lmerTest)   # For p-values in mixed models
library(emmeans)    # For post-hoc tests
library(akima)      # For interpolation in heatmaps
library(viridis)    # For color scales

# Import Data -------------------------------------------------------------
# Define repository path - update this to your path when sharing
DATA_PATH <- "PATH_TO_YOUR_REPOSITORY/data"

# Load previously identified significant sites
sig_sites_data <- read.csv(file.path(DATA_PATH, 'significant_mua_sites/targacq_-1to1_sig01.csv'))
significant_sites <- unique(sig_sites_data$sites)

# CCA Analysis across delays and timepoints -------------------------------
# Optional: If you've already run the analysis and saved the results,
# you can skip directly to the "Load previously run results" section

# Parameters for CCA analysis
drugs <- c('oxytocin', 'saline')
days_data <- sig_sites_data %>% filter(drugs != 'unspecified') %>% select(days) %>% distinct(days)
days <- days_data$days
outcomes <- c('other', 'none')
trialtypes <- c('choice')
states <- c('high', 'low')
administrations <- c('pre', 'post')
time_points <- seq(-500, 300, by = 25)  # base time points in ms, adjust as needed
time_delays <- seq(-100, 100, by = 5)   # delays in ms, adjust as needed

# Create list to store all results
all_cca_results <- list()

# Perform CCA analysis between BLA and ACC across different time delays
for (d in 1:length(days)){
  cur_day <- days[d]
  print(paste("Processing day:", d, "of", length(days)))
  
  for (t in 1:length(outcomes)){
    cur_outcome <- outcomes[t]
    
    for (a in 1:length(administrations)) {
      cur_admin <- administrations[a]
      
      # First normalize across all data
      normalized_data <- early_sites_psth %>%
        mutate(state = ifelse(days %in% high_combo_prefindex_days, 'high', 'low')) %>%
        filter(drugs != 'unspecified', outcomes != 'errors' & outcomes != 'self' & outcomes != 'both',
               days == cur_day, administration == cur_admin, outcomes == cur_outcome) %>%
        unite(day_block_trial, c(days, blocks, trials), remove = F) %>%
        group_by(channel, region) %>%
        mutate(across(matches("^[0-9]"), ~scale(.))) %>%
        ungroup()
      
      # Skip if no data for this condition
      if(nrow(normalized_data) == 0) {
        print(paste("No data for day:", cur_day, "outcome:", cur_outcome, 
                    "admin:", cur_admin))
        next
      }
      
      # Get metadata for this condition
      metadata <- normalized_data %>%
        select(days, drugs, state, administration, outcomes, trialtypes, monkeys) %>%
        distinct()
      
      # Loop through base time points
      for(base_time in time_points) {
        # Loop through delays
        for(delay in time_delays) {
          # Convert delay from ms to column value
          delay_s <- delay
          
          if(delay >= 0) {
            bla_time <- base_time
            acc_time <- base_time + delay_s
          } else {
            acc_time <- base_time
            bla_time <- base_time - delay_s
          }
          
          # Convert times to column names
          bla_col <- as.character(bla_time)
          acc_col <- as.character(acc_time)
          
          # Check if these time points exist in the data
          if(!(bla_col %in% names(normalized_data)) || !(acc_col %in% names(normalized_data))) {
            next
          }
          
          tryCatch({
            # Create BLA matrix
            bla_matrix <- normalized_data %>%
              filter(region == "bla") %>%
              select(day_block_trial, channel, all_of(bla_col)) %>%
              pivot_wider(names_from = channel,
                          values_from = all_of(bla_col)) %>%
              select(-day_block_trial) %>%
              as.matrix()
            
            # Create ACC matrix
            acc_matrix <- normalized_data %>%
              filter(region == "acc") %>%
              select(day_block_trial, channel, all_of(acc_col)) %>%
              pivot_wider(names_from = channel,
                          values_from = all_of(acc_col)) %>%
              select(-day_block_trial) %>%
              as.matrix()
            
            # Run CCA
            cca_result <- cancor(bla_matrix, acc_matrix)
            
            # Store result with all metadata
            result <- list(
              correlations = cca_result$cor,
              bla_coeffs = cca_result$xcoef,
              acc_coeffs = cca_result$ycoef,
              bla_channels = colnames(bla_matrix),
              acc_channels = colnames(acc_matrix),
              n_trials = nrow(bla_matrix),
              metadata = metadata,
              condition = list(
                day = cur_day,
                drug = metadata$drugs[1],
                state = metadata$state[1],
                administration = cur_admin,
                outcome = cur_outcome,
                trialtype = metadata$trialtypes[1],
                monkey = metadata$monkeys[1],
                delay = delay,
                base_time = base_time,
                bla_time = bla_time,
                acc_time = acc_time
              )
            )
            
            # Add to results list
            all_cca_results[[length(all_cca_results) + 1]] <- result
            
          }, error = function(e) {
            print(paste("Error in CCA for day:", cur_day, 
                        "outcome:", cur_outcome, 
                        "admin:", cur_admin,
                        "delay:", delay, "ms",
                        "base time:", base_time))
          })
        }
      }
    }
  }
}

# Create summary data frame with delay information
summary_df <- do.call(rbind, lapply(all_cca_results, function(result) {
  data.frame(
    day = result$condition$day,
    drug = result$condition$drug,
    state = result$condition$state,
    administration = result$condition$administration,
    outcome = result$condition$outcome,
    trialtype = result$condition$trialtype,
    monkey = result$condition$monkey,
    delay = result$condition$delay,
    bla_time = result$condition$bla_time,
    acc_time = result$condition$acc_time,
    first_canonical_corr = result$correlations[1],
    n_trials = result$n_trials,
    n_bla_channels = length(result$bla_channels),
    n_acc_channels = length(result$acc_channels),
    base_time = result$condition$base_time
  )
}))

# Save results
write.csv(summary_df, file.path(DATA_PATH, 'mua/CCA/cca_results_summary.csv'))

# Load previously run results ---------------------------------------------
# If you've already run the CCA analysis and saved results
summary_df <- read.csv(file.path(DATA_PATH, 'mua/CCA/cca_results_summary.csv'))

# Post - Pre CCA Differences ----------------------------------------------
# Calculate differences in canonical correlations between post and pre drug administration
cca_difference <- summary_df %>%
  # First group by day to get post-pre difference for each day
  group_by(drug, state, outcome, base_time, delay, day) %>%
  summarise(
    diff_correlation = mean(first_canonical_corr[administration == "post"], na.rm = TRUE) - 
      mean(first_canonical_corr[administration == "pre"], na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Then get the mean and SEM across days
  group_by(drug, state, outcome, base_time, delay) %>%
  summarise(
    mean_diff = mean(diff_correlation, na.rm = TRUE),
    sem_diff = sd(diff_correlation, na.rm = TRUE)/sqrt(n()),
    .groups = 'drop'
  ) %>%
  ungroup() %>%
  # Rename to match visualization code
  rename(diff_correlation = mean_diff) %>%
  # Filter to focus on relevant delays
  filter(delay >= -50 & delay <= 50)

# Create interpolation function for smoother heatmaps
interpolate_facet_data <- function(data) {
  # Create a regular grid for interpolation
  x_seq <- seq(min(data$base_time), max(data$base_time), length.out = 50)
  y_seq <- seq(min(data$delay), max(data$delay), length.out = 50)
  
  # Perform bilinear interpolation
  interp_result <- interp(
    x = data$base_time,
    y = data$delay,
    z = data$mean_correlation,
    xo = x_seq,
    yo = y_seq,
    linear = TRUE,  # Use smooth interpolation
    extrap = FALSE   # Don't extrapolate beyond data bounds
  )
  
  # Convert interpolated results to data frame
  expand.grid(base_time = interp_result$x,
              delay = interp_result$y) %>%
    mutate(mean_correlation = as.vector(interp_result$z))
}

# Interpolate each facet separately for smoother visualization
interpolated_data <- cca_difference %>%
  group_by(drug, state, outcome) %>%
  group_modify(~interpolate_facet_data(data.frame(
    base_time = .x$base_time,
    delay = .x$delay,
    mean_correlation = .x$diff_correlation
  ))) %>%
  ungroup()

# Create heatmap visualization
cca_heatmap <- ggplot(interpolated_data, 
                      aes(x = base_time, y = delay, fill = mean_correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0,
                       limits = c(-0.2, 0.2), oob = scales::squish,
                       name = "Change in\nCanonical\nCorrelation\n(Post - Pre)") +
  facet_grid(state ~ drug + outcome) +
  labs(x = "Time from Target Acquisition (ms)", 
       y = "Delay (ms)\nNegative = ACC leads, Positive = BLA leads",
       title = "Change in BLA-ACC Correlation After Drug Administration") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 1) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "grey", alpha = 1) +
  geom_hline(yintercept = -15, linetype = "dashed", color = "grey", alpha = 1) +
  geom_vline(xintercept = -160, linetype = 'dashed', color = 'black', alpha = 1) +
  geom_vline(xintercept = 10, linetype = 'dashed', color = 'black', alpha = 1) +
  xlim(-500, 300)

print(cca_heatmap)

# Time course of changes with error ribbons
# Separate analysis by direction (ACC→BLA vs BLA→ACC)
plot_data <- cca_difference %>% 
  filter(base_time > -500 & base_time < 300, 
         (delay > 10 & delay < 60) | (delay < -10 & delay > -60)) %>%
  mutate(direction = ifelse(delay < 0, 'ACC leads BLA', 'BLA leads ACC'))

# Create line plot with error ribbons
cca_time_plot <- ggplot(plot_data %>%
                          group_by(drug, state, outcome, base_time, direction) %>%
                          summarise(
                            mean_change = mean(diff_correlation, na.rm = TRUE),
                            se_change = sem_diff,
                            .groups = 'drop'
                          ),
                        aes(x = base_time, y = mean_change, 
                            color = interaction(outcome, drug),
                            fill = interaction(outcome, drug))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_change - se_change, 
                  ymax = mean_change + se_change),
              alpha = 0.2, linetype = 0) +
  facet_wrap(~state + direction, ncol = 2, scales = 'free') +
  scale_color_manual(values = c(
    "none.oxytocin" = "#D4787D",      # Rose
    "other.oxytocin" = "#4C9085",     # Teal
    "none.saline" = "#E5B0B5",        # Light rose
    "other.saline" = "#93B5B3"        # Light teal
  )) +
  scale_fill_manual(values = c(
    "none.oxytocin" = "#D4787D",      # Rose
    "other.oxytocin" = "#4C9085",     # Teal
    "none.saline" = "#E5B0B5",        # Light rose
    "other.saline" = "#93B5B3"        # Light teal
  )) +
  theme_classic() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = .5) +
  geom_vline(xintercept = -150, linetype = 'dashed', color = 'black', alpha = .5) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', alpha = .5) +
  labs(x = "Time from Target Acquisition (ms)",
       y = "Change in Correlation",
       color = "Drug/Outcome",
       fill = "Drug/Outcome") +
  theme(legend.position = "bottom") +
  ylim(-.1,.2)

print(cca_time_plot)

# Statistical analysis: Mixed effects model
lme_data <- cca_difference %>% 
  filter(base_time > -500 & base_time < 300, 
         (delay > 10 & delay < 60) | (delay < -10 & delay > -60)) %>%
  mutate(direction = ifelse(delay < 0, 'ACC leads BLA', 'BLA leads ACC')) %>%
  group_by(drug, state, monkey, day, outcome, base_time, direction) %>%
  summarise(mean_cca_diff = mean(diff_correlation, na.rm = TRUE), 
            .groups = 'drop') %>%
  mutate(drug = factor(drug, levels = c("oxytocin", "saline"))) %>%
  mutate(state = factor(state, levels = c("high", "low"))) %>%
  mutate(direction = factor(direction, levels = c('BLA leads ACC', 'ACC leads BLA')))

# Run the mixed-effects model
mixed_model_cca <- lmer(mean_cca_diff ~ drug * state * direction * outcome + 
                          (1|day) + (1|monkey), 
                        data = lme_data)

# Summary of the mixed model
summary(mixed_model_cca)

# Get estimated marginal means and comparisons
emmeans_comb <- emmeans(mixed_model_cca, ~ drug * state * outcome * direction)
pairs(emmeans_comb, adjust = "tukey")

# Direction Dominance Index (DDI) Analysis -------------------------------
# Calculate DDI to measure directional dominance in BLA-ACC communication

# Define DDI calculation function
calculate_ddi <- function(data) {
  # Calculate DDI for each time point
  ddi_by_time <- data %>%
    group_by(base_time) %>%
    summarise(
      # Calculate mean correlations for each direction at each time point
      CFF_t = mean(first_canonical_corr[delay > 0], na.rm = TRUE),
      CFB_t = mean(first_canonical_corr[delay < 0], na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Calculate normalization term components
  N_t <- n_distinct(data$base_time)
  N_dt <- n_distinct(data$delay)
  C_sum <- sum(data$first_canonical_corr, na.rm = TRUE)
  norm_term <- C_sum/sqrt(N_t * N_dt)
  
  # Calculate final DDI
  ddi_by_time %>%
    mutate(
      DDI = (CFF_t - CFB_t)/norm_term
    )
}

# Define analysis window (10-60ms absolute delay)
min_delay <- 10
max_delay <- 60

# Calculate DDI for each day and administration separately
ddi_by_day_admin <- summary_df %>%
  filter((delay > min_delay & delay < max_delay) | 
           (delay < -min_delay & delay > -max_delay)) %>%
  group_by(drug, state, outcome, day, administration) %>%
  group_modify(~calculate_ddi(.x)) %>%
  ungroup()

# Calculate post-pre difference in DDI
ddi_difference <- ddi_by_day_admin %>%
  group_by(drug, state, outcome, day, base_time) %>%
  summarise(
    diff_DDI = DDI[administration == "post"] - DDI[administration == "pre"],
    .groups = 'drop'
  ) %>%
  ungroup()

# Calculate mean and SEM across days
ddi_summary <- ddi_difference %>%
  group_by(drug, state, outcome, base_time) %>%
  summarise(
    mean_DDI = mean(diff_DDI, na.rm = TRUE),
    sem_DDI = sd(diff_DDI, na.rm = TRUE)/sqrt(n()),
    .groups = 'drop'
  ) %>%
  ungroup()

# Plot DDI time course with error ribbons
ddi_plot <- ggplot(ddi_summary, aes(x = base_time, y = mean_DDI, 
                                    color = interaction(drug, state),
                                    fill = interaction(drug, state))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = -150, linetype = 'dashed', color = 'grey', alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey', alpha = 0.5) +
  geom_ribbon(aes(ymin = mean_DDI - sem_DDI,
                  ymax = mean_DDI + sem_DDI), linetype = 0, alpha = 0.25) +
  geom_line(size = .8, alpha = .9) +
  scale_color_manual(values = c(
    "oxytocin.high" = "#ec8249",
    "saline.high" = "#1e71a8",
    "oxytocin.low" = "#EDB120",
    "saline.low" = "#4DBEEE"
  )) +
  scale_fill_manual(values = c(
    "oxytocin.high" = "#ec8249",
    "saline.high" = "#1e71a8",
    "oxytocin.low" = "#EDB120",
    "saline.low" = "#4DBEEE"
  )) +
  facet_wrap(~state + outcome, ncol = 1, scales = 'free') +
  labs(x = "Time from Target Acquisition (ms)",
       y = "Direction Dominance Index (Post-Pre)",
       title = "Change in Direction Dominance Index Over Time",
       subtitle = "Positive = increased BLA→ACC dominance, Negative = increased ACC→BLA dominance") +
  theme_classic() +
  xlim(-500, 300)

print(ddi_plot)

# Analysis of DDI during choice formation period -------------------------
# First get the window DDI values during choice period (-150 to 0 ms)
ddi_window <- ddi_difference %>%
  filter(base_time >= -150 & base_time <= 0) %>%
  group_by(drug, state, outcome, day) %>%
  summarise(mean_DDI = mean(diff_DDI, na.rm = TRUE), 
            .groups = 'drop') %>%
  ungroup()

# Create box plot for choice period
ddi_boxplot <- ggplot(ddi_window, aes(x = interaction(outcome), y = mean_DDI, 
                                      fill = interaction(drug, state))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(alpha = .5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             size = 1.5, alpha = 0.8, aes(color = interaction(drug, state))) +
  scale_fill_manual(values = c(
    "oxytocin.high" = "#ec8249",
    "saline.high" = "#1e71a8",
    "oxytocin.low" = "#EDB120",
    "saline.low" = "#4DBEEE"
  )) +
  scale_color_manual(values = c(
    "oxytocin.high" = "#ec8249",
    "saline.high" = "#1e71a8",
    "oxytocin.low" = "#EDB120",
    "saline.low" = "#4DBEEE"
  )) +
  facet_wrap(~state, ncol = 2, scales = 'free') +
  labs(x = "Outcome",
       y = "Change in Direction Dominance Index (Post-Pre)",
       title = "DDI Change During Choice Formation (-150 to 0 ms)",
       subtitle = "Positive = increased BLA→ACC dominance, Negative = increased ACC→BLA dominance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(ddi_boxplot)

# Calculate summary statistics for bar plot
ddi_summary <- ddi_window %>%
  group_by(outcome, drug, state) %>%
  summarise(
    mean_ddi = mean(mean_DDI, na.rm = TRUE),
    sem = sd(mean_DDI, na.rm = TRUE)/sqrt(n()),
    .groups = 'drop'
  )

# Create bar plot for clearer visualization
ddi_summary$outcome <- factor(ddi_summary$outcome, levels = c('other', 'none'))
ddi_barplot <- ggplot(ddi_summary, 
                      aes(x = outcome, y = mean_ddi, 
                          fill = interaction(drug, state))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.9),
           color = "black",
           width = 0.8) +
  geom_errorbar(aes(ymin = mean_ddi - sem,  
                    ymax = mean_ddi + sem),  
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(
    "oxytocin.high" = "#ec8249",
    "saline.high" = "#1e71a8",
    "oxytocin.low" = "#EDB120",
    "saline.low" = "#4DBEEE"
  )) +
  labs(x = "Outcome",
       y = "Change in Direction Dominance Index (Post-Pre)",
       title = "DDI Change During Choice Formation (-150 to 0 ms)",
       subtitle = "Positive = increased BLA→ACC dominance, Negative = increased ACC→BLA dominance",
       fill = "Drug/State") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(ddi_barplot)

# Statistical analysis of DDI differences
# Wilcoxon tests to compare each condition to zero and between drugs
stats_summary <- ddi_window %>%
  # First get stats for each drug/state/outcome combo
  group_by(state, outcome, drug) %>%
  summarise(
    # Test if different from zero
    zero_W_stat = tryCatch(
      wilcox.test(mean_DDI, mu = 0, exact = FALSE)$statistic,
      error = function(e) NA
    ),
    zero_p_value = tryCatch(
      wilcox.test(mean_DDI, mu = 0, exact = FALSE)$p.value,
      error = function(e) NA
    ),
    mean_ddi = mean(mean_DDI, na.rm = TRUE),
    std_ddi = sd(mean_DDI, na.rm = TRUE),
    n = n(),
    .groups = 'keep'
  ) %>%
  # Pivot wider to compare drugs
  pivot_wider(
    id_cols = c(state, outcome),
    names_from = drug,
    values_from = c(zero_W_stat, zero_p_value, mean_ddi, std_ddi, n)
  ) %>%
  # Add drug comparison
  rowwise() %>%
  mutate(
    # Get the drug comparison test
    drug_W_stat = tryCatch(
      wilcox.test(
        ddi_window$mean_DDI[ddi_window$drug == "oxytocin" & 
                              ddi_window$state == state & 
                              ddi_window$outcome == outcome],
        ddi_window$mean_DDI[ddi_window$drug == "saline" & 
                              ddi_window$state == state & 
                              ddi_window$outcome == outcome],
        exact = FALSE
      )$statistic,
      error = function(e) NA
    ),
    drug_p_value = tryCatch(
      wilcox.test(
        ddi_window$mean_DDI[ddi_window$drug == "oxytocin" & 
                              ddi_window$state == state & 
                              ddi_window$outcome == outcome],
        ddi_window$mean_DDI[ddi_window$drug == "saline" & 
                              ddi_window$state == state & 
                              ddi_window$outcome == outcome],
        exact = FALSE
      )$p.value,
      error = function(e) NA
    )
  ) %>%
  ungroup()
