# OT-BLA-ACC

## Repository Structure

The repository consists of three main analysis scripts:

1. **Behavioral Analysis** (`behavior_analyses.R`): Analyzes choice behavior, reaction time, and incomplete trial rates during social reward allocation task.
2. **Neural Firing Rate Analysis** (`firing_rate_analysis.R`): Examines how oxytocin affects neural activity in BLA and ACC regions.
3. **Inter-regional Communication Analysis** (`cca_analysis.R`): Implements canonical correlation analysis to quantify directional information flow between brain regions.

## Script Descriptions

### 1. Behavioral Analysis (`behavioral_analysis.R`)
Analyzes choice behavior in a prosocial decision-making task where subjects chose between different reward outcomes (self, both, other, none).

**Key Analyses:**
- Preference index to capture prosocial choice rates
- Incomplete trial (a.k.a. Error) rate analysis across session
- Reaction time analysis
- Linear mixed effects models of each


### 2. Neural Analysis (`neural_analysis.R`)
Analyzes neural firing rates in BLA and ACC regions before and after oxytocin administration.

**Key Analyses:**
- Neural activity time-locked to target acquisition
- Mixed effects models comparing drug conditions
- Time bin-specific comparisons

### 3. CCA Analysis (`cca_analysis.R`)
Implements canonical correlation analysis to measure information flow between BLA and ACC regions across different time delays.

**Key Analyses:**
- Canonical correlation computation across time delays
- Post-Pre oxytocin effect visualization
- Direction Dominance Index (DDI) calculation
- Statistical analysis of directional information flow

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
