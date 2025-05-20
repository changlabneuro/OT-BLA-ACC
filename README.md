# OT-BLA-ACC

## Repository Structure

The repository consists of three main analysis scripts:

1. **Behavioral Analysis** (`behavior_analyses.R`): Analyzes choice behavior, reaction time, and incomplete trial rates during social reward allocation task.
2. **Neural Firing Rate Analysis** (`firing_rate_analysis.R`): Examines how oxytocin affects neural activity in BLA and ACC regions.
3. **Inter-regional Communication Analysis** (`cca_analysis.R`): Implements canonical correlation analysis to quantify directional information flow between brain regions.

## Data Requirements

The scripts require the following data files, organized within a `data` directory:

### Behavioral Data
- `behavioral_data.csv`: Main behavioral dataset (excluding error trials)
- `behavioral_data_with_errors.csv`: Complete behavioral dataset including error trials
- `high_prosocial_days.csv`: List of experimental days classified as high prosocial baseline days

### Neural Activity Data
- `targacq_mua.csv`: Multi-unit activity data time-locked to target acquisition
- `significant_mua_sites/targacq_-1to1_sig01.csv`: List of recording sites showing significant activity modulation
- `mua/normalized_to_cueon/targacq_nas_corrected/*.csv`: Normalized MUA data files by day

### CCA Analysis Data
- `mua/directionality_timebins/targacq_raw/binned_psth_5ms_-1to1.csv`: Binned neural activity for directional analysis
- `mua/directionality_timebins/targacq_raw/raw_psth_labels.csv`: Labels for binned neural data

## Script Descriptions

### 1. Behavioral Analysis (`behavioral_analysis.R`)
Analyzes choice behavior in a prosocial decision-making task where subjects chose between different reward outcomes (self, both, other, none).

**Key Analyses:**
- Preference index to capture prosocial choice rates
- Incomplete trial (a.k.a. Error) rate analysis across session
- Reaction time analysis
- Linear mixed effects models of each

**Data Requirements:** `behavioral_data.csv`, `behavioral_data_with_errors.csv`, `high_prosocial_days.csv`

### 2. Neural Analysis (`neural_analysis.R`)
Analyzes neural firing rates in BLA and ACC regions before and after oxytocin administration.

**Key Analyses:**
- Neural activity time-locked to target acquisition
- Mixed effects models comparing drug conditions
- Time bin-specific comparisons

**Data Requirements:** `targacq_mua.csv`, `significant_mua_sites/targacq_-1to1_sig01.csv`

### 3. CCA Analysis (`cca_analysis.R`)
Implements canonical correlation analysis to measure information flow between BLA and ACC regions across different time delays.

**Key Analyses:**
- Canonical correlation computation across time delays
- Post-Pre oxytocin effect visualization
- Direction Dominance Index (DDI) calculation
- Statistical analysis of directional information flow

**Data Requirements:** `high_prosocial_days.csv`, `mua/directionality_timebins/targacq_raw/binned_psth_5ms_-1to1.csv`, `mua/directionality_timebins/targacq_raw/raw_psth_labels.csv`

## Dependencies

The scripts require the following R packages:

```r
# Main packages
library(tidyverse)  # For data manipulation and visualization
library(lme4)       # For linear mixed effects models
library(lmerTest)   # For p-values in mixed models
library(emmeans)    # For post-hoc tests

# Additional packages
library(akima)      # For interpolation (CCA analysis)
library(viridis)    # For color scales
```

## Citation

If you use this code in your research, please cite:

[Citation information will be added upon publication]
