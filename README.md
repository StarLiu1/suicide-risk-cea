# Suicide Risk Prediction Cost-Effectiveness Model

Replication of the Ross et al. (2021) economic evaluation model for suicide risk prediction using R and the `hesim` package.

## Overview

This project implements a discrete-time state-transition Markov model to evaluate the cost-effectiveness of suicide risk prediction interventions in US primary care patients. The model compares three strategies:

- **No Prediction** - Usual care (baseline)
- **ACF Intervention** - Active Contact and Follow-up (RR = 0.83, cost = $96/year)
- **CBT Intervention** - Cognitive Behavioral Therapy (RR = 0.47, cost = $1,088/year)

## Model Structure

The model uses a 3-state structure replicated across 1000 risk strata:

<img width="701" height="316" alt="image" src="https://github.com/user-attachments/assets/5dde7b02-f51b-4305-a830-3dc293496aaa" />


### Health States

| State | Description |
|-------|-------------|
| **No_Attempt (1)** | Patient has never made a suicide attempt |
| **Prior_Attempt (2)** | Patient has survived at least one suicide attempt |
| **Dead (3)** | Absorbing state (suicide or non-suicide mortality) |

### Key Transitions

- **No_Attempt → Prior_Attempt**: Suicide attempt survived
- **No_Attempt → Dead**: Suicide death or non-suicide mortality
- **Prior_Attempt → Prior_Attempt**: Subsequent attempt (survived) or no new attempt
- **Prior_Attempt → Dead**: Suicide death or non-suicide mortality

## Target Rates (from Ross et al. 2021)

| Parameter | Value |
|-----------|-------|
| Suicide attempts | 175 per 100,000 person-years |
| Suicide deaths | 15 per 100,000 person-years |
| Death per attempt | 8.81% |
| Prior attempt multiplier | 1.54× |
| Base utility (EQ-5D) | 0.866 |
| Mean population age | 48.8 years (SD 17.2) |

## Project Structure

```
suicide-risk-cea/
├── R/
│   ├── 00-validation-checks.R    # Transition matrix validation tests
│   ├── 00-diagnostics.R          # Simulation diagnostic functions
│   ├── 01-model-setup.R          # hesim data structures, risk strata, patients
│   ├── 02-parameters.R           # Model parameters, define_model(), tparams
│   ├── 03-transitions.R          # CohortDtstmTrans transition model
│   ├── 04-costs-utilities.R      # StateVals for costs and utilities
│   ├── 05-simulation.R           # Main simulation engine (hybrid approach)
│   ├── 05-simulation2.R          # Alternative simulation implementation
│   ├── 06-risk-prediction.R      # Sensitivity/specificity targeting layer
│   ├── 07-complete-lifetime-cea.R           # Full lifetime CEA simulation
│   └── 07-complete-lifetime-cea-vectorized.R # Optimized vectorized version
├── External Data/
│   ├── lifetable.xlsx            # CDC life table mortality data
│   └── age_specific_suicide.xlsx # Age-specific suicide rates
├── data/                         # Generated .RData files
│   ├── hesim_setup.RData
│   ├── hesim_parameters.RData
│   ├── hesim_transitions.RData
│   └── hesim_costs_utilities.RData
├── output/
│   ├── results/                  # Simulation results (.RData)
│   └── figures/                  # Generated plots (.png)
├── docs/
│   └── state_transition_diagram.png
└── README.md
```

## Key Features

- **Individual patient tracking** with 25,000-100,000 simulated patients
- **1000 risk strata** calibrated to match paper's logit-normal distribution
- **Age-dependent mortality** using CDC life tables
- **Age-stratified costs** (18-44, 45-64, 65+ years)
- **Risk prediction targeting** with configurable sensitivity/specificity
- **Vectorized simulation** for performance (~2 seconds per strategy)
- **Lifetime horizon** (60+ cycles with discounting at 3%)

## Requirements

```r
# Core packages
install.packages(c("hesim", "data.table", "ggplot2", "readxl"))

# For development
install.packages(c("devtools", "usethis"))
```

- R ≥ 4.0
- `hesim` - Health economic simulation and CEA
- `data.table` - Efficient data manipulation
- `ggplot2` - Visualization
- `readxl` - Reading external mortality data

## Usage

```r
# Run scripts in order:
source("R/01-model-setup.R")
source("R/02-parameters.R")
source("R/03-transitions.R")
source("R/04-costs-utilities.R")
source("R/05-simulation.R")

# For risk prediction analysis:
source("R/06-risk-prediction.R")

# For complete lifetime CEA:
source("R/07-complete-lifetime-cea-vectorized.R")
```

## Risk Prediction Parameters

The model implements risk prediction with configurable accuracy:

| Parameter | Base Case (Table 2) |
|-----------|---------------------|
| Specificity | 95% |
| Sensitivity | 25% |
| PPV (attempts) | ~0.8% (ACF), ~1.7% (CBT) |

## Cost-Effectiveness Results

At $150,000/QALY threshold with 95% specificity:
- **ACF**: Cost-effective at ≥17% sensitivity
- **CBT**: Cost-effective at ≥36% sensitivity

## Reference

Ross EL, Zuromski KL, Reis BY, Nock MK, Kessler RC, Smoller JW. Accuracy Requirements for Cost-effective Suicide Risk Prediction Among Primary Care Patients in the US. *JAMA Psychiatry*. 2021;78(6):642-650. doi:[10.1001/jamapsychiatry.2021.0089](https://doi.org/10.1001/jamapsychiatry.2021.0089)

## License

MIT
