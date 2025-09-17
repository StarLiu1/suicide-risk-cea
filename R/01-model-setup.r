# 1. Hesim Setup and Data Structure
# File: R/01-hesim-setup.R

library(hesim)
library(data.table)

# Clear workspace
rm(list = ls())

cat("Setting up hesim data structures...\n")

# Define strategies (interventions)
strategies <- data.table(
  strategy_id = 1:3,
  strategy_name = c("No_Prediction", "ACF_Intervention", "CBT_Intervention")
)

cat("Strategies defined:\n")
print(strategies)

# Define patients (simplified cohort)
# For Phase 1, we'll use a single representative patient per risk stratum
n_patients <- 1  # One representative patient per risk stratum in Phase 1
patients <- data.table(
  patient_id = 1:n_patients,
  age = 48.8,  # Mean age from paper
  sex = "Both"  # Mixed population
)

cat("\nPatients defined:\n")
print(patients)

# Define health states
# We need states for each risk stratum: No attempts, Prior attempt, Dead
n_risk_strata <- 10

# Create states: each risk stratum has 3 health states
states <- data.table()
state_id <- 1

for (stratum in 1:n_risk_strata) {
  # Add the 3 states for this stratum
  stratum_states <- data.table(
    state_id = state_id:(state_id + 2),
    stratum = stratum,
    state_type = c("no_attempts", "prior_attempt", "dead"),
    state_name = paste0("s", stratum, "_", c("no_attempts", "prior_attempt", "dead"))
  )
  
  states <- rbind(states, stratum_states)
  state_id <- state_id + 3
}

cat("\nStates defined (showing first 12):\n")
print(head(states, 12))
cat("Total states:", nrow(states), "\n")

# Define risk stratum characteristics
# This will be used to parameterize transition probabilities
risk_strata <- data.table(
  stratum_id = 1:n_risk_strata,
  # Use calibrated rates from Ross et al. paper
  baseline_attempt_rate = c(
    rep(0.00008, 6),     # 0-60th percentile (very low risk)
    rep(0.001, 2),       # 60-80th percentile (low risk) 
    rep(0.005, 1),       # 80-90th percentile (medium risk)
    rep(0.15437, 1)      # >90th percentile (high risk)
  ),
  population_weight = c(
    rep(0.10, 6),        # 60% of population in very low risk
    rep(0.10, 2),        # 20% of population in low risk
    rep(0.10, 1),        # 10% of population in medium risk  
    rep(0.10, 1)         # 10% of population in high risk
  )
)

# Verify population weights sum to 1
total_weight <- sum(risk_strata$population_weight)
cat("\nRisk strata population weights sum to:", total_weight, "\n")

if (abs(total_weight - 1.0) > 1e-6) {
  warning("Population weights don't sum to 1.0!")
}

cat("Risk strata defined:\n")
print(risk_strata)

# Create hesim_data object
# This combines all the structural elements
hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients, 
  states = states
)

cat("\nHesim data object created successfully!\n")
cat("- Strategies:", nrow(hesim_dat$strategies), "\n")
cat("- Patients:", nrow(hesim_dat$patients), "\n") 
cat("- States:", nrow(hesim_dat$states), "\n")

# Create expanded dataset for modeling
# This creates all combinations of strategies, patients, and states
expanded_dat <- expand(hesim_dat)

cat("\nExpanded dataset created:\n")
cat("- Total rows:", nrow(expanded_dat), "\n")
cat("- Combinations: strategies x patients x states =", 
    nrow(strategies), "x", nrow(patients), "x", nrow(states), "=",
    nrow(strategies) * nrow(patients) * nrow(states), "\n")

# Save objects for other scripts
save(
  strategies, patients, states, risk_strata, hesim_dat, expanded_dat,
  n_risk_strata, n_patients,
  file = "data/hesim_setup.RData"
)

cat("\n✓ Setup complete! Data saved to data/hesim_setup.RData\n")
cat("\nNext: Run 02-hesim-parameters.R to define model parameters\n")