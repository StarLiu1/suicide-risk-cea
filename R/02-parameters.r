# 2. Hesim Model Parameters  
# File: R/02-hesim-parameters.R

library(hesim)
library(data.table)

# Load setup data
load("data/hesim_setup.RData")

cat("Defining model parameters for hesim...\n")

# Model timing parameters
n_cycles <- 80      # Reduced for Phase 1 (vs 80 for lifetime in full model)
cycle_length <- 1   # 1 year cycles
discount_rate <- 0.03  # 3% annual discount rate

cat("Model timing:\n")
cat("- Cycles:", n_cycles, "\n")
cat("- Cycle length:", cycle_length, "year\n")
cat("- Discount rate:", discount_rate, "\n")

# Intervention parameters
intervention_params <- list(
  # Relative risk of suicide attempt (from meta-analyses in paper)
  rr = c(
    "No_Prediction" = 1.0,      # Baseline
    "ACF_Intervention" = 0.83,   # Active Contact and Follow-up
    "CBT_Intervention" = 0.47    # Cognitive Behavioral Therapy
  ),
  
  # Annual intervention costs (2016 USD)
  annual_cost = c(
    "No_Prediction" = 0,
    "ACF_Intervention" = 96,     # Safety planning + telephone follow-up
    "CBT_Intervention" = 1088    # Individual CBT sessions
  ),
  
  # Intervention uptake rates
  uptake = c(
    "No_Prediction" = 1.0,       # No intervention, so 100% "uptake"
    "ACF_Intervention" = 0.994,  # 99.4% uptake from Stanley et al.
    "CBT_Intervention" = 0.899   # 89.9% uptake from Rudd et al.
  )
)

cat("\nIntervention parameters defined:\n")
print(data.frame(
  Strategy = names(intervention_params$rr),
  Relative_Risk = intervention_params$rr,
  Annual_Cost = intervention_params$annual_cost,
  Uptake_Rate = intervention_params$uptake
))

# Clinical parameters
clinical_params <- list(
  # Death probability per suicide attempt (from paper)
  death_per_attempt = 0.0881,  # 8.81%
  
  # Prior attempt effect (54% of attempts are from those with prior attempts)
  prior_attempt_multiplier = 1.54,  # 54% higher rate for those with prior attempts
  
  # Base utility (EQ-5D for primary care population)
  base_utility = 0.866
)

cat("\nClinical parameters:\n")
cat("- Death per attempt:", clinical_params$death_per_attempt, "\n")
cat("- Prior attempt multiplier:", clinical_params$prior_attempt_multiplier, "\n") 
cat("- Base utility:", clinical_params$base_utility, "\n")

# Cost parameters (2016 USD)
cost_params <- list(
  # Suicide attempt costs
  nonfatal_attempt_medical = 10830,    # Medical cost per nonfatal attempt
  nonfatal_attempt_productivity = 17369, # Productivity cost per nonfatal attempt
  fatal_attempt_medical = 4354,        # Medical cost per fatal attempt (age-adjusted)
  fatal_attempt_productivity = 61150,  # Productivity cost per fatal attempt (age-adjusted)
  
  # Background healthcare costs by age group (annual)
  bg_medical_18_44 = 4016,
  bg_medical_45_64 = 7648,
  bg_medical_65plus = 11740,
  
  # Risk assessment cost (evaluation after positive screen)
  evaluation_cost = 76
)

cat("\nCost parameters (2016 USD):\n")
cat("- Nonfatal attempt medical:", cost_params$nonfatal_attempt_medical, "\n")
cat("- Fatal attempt medical:", cost_params$fatal_attempt_medical, "\n")
cat("- Background medical (45-64):", cost_params$bg_medical_45_64, "\n")

# Function to create transition probability parameters for hesim
create_transition_params <- function() {
  
  cat("\nCreating transition probability parameters...\n")
  
  # We need to create parameters for each state transition
  # For hesim, we'll use logit models for transition probabilities
  
  # Create parameter table for all transitions
  n_states <- nrow(states)
  n_strategies <- nrow(strategies)
  n_pts <- nrow(patients)
  
  # Create expanded dataset for transitions
  trans_data <- expand(hesim_dat, by = c("strategies", "patients"))
  
  # Add stratum information based on patient_id 
  # (in Phase 1, patient_id corresponds to risk stratum)
  trans_data[, stratum := ((patient_id - 1) %% n_risk_strata) + 1]
  
  # Add baseline attempt rates
  trans_data <- merge(trans_data, risk_strata[, .(stratum_id, baseline_attempt_rate)], 
                      by.x = "stratum", by.y = "stratum_id")
  
  # Add intervention effects
  trans_data[, rr := intervention_params$rr[strategy_name]]
  
  # Calculate adjusted attempt rates
  trans_data[, adjusted_attempt_rate := baseline_attempt_rate * rr]
  
  cat("Transition parameters dataset created with", nrow(trans_data), "rows\n")
  
  return(trans_data)
}

# Create the transition parameters
transition_data <- create_transition_params()

cat("\nSample of transition data:\n")
print(head(transition_data))

# For hesim, we need to define the transition model structure
# We'll create this in the next script, but prepare the parameters here

# Save all parameters
save(
  n_cycles, cycle_length, discount_rate,
  intervention_params, clinical_params, cost_params,
  transition_data,
  file = "data/hesim_parameters.RData"
)

cat("\n✓ Parameters defined and saved to data/hesim_parameters.RData\n")
cat("\nNext: Run 03-hesim-transitions.R to create transition model\n")