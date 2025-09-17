# 2. Hesim Model Parameters - PHASE 2 INDIVIDUAL PATIENT MODEL
# File: R/02-parameters.R
# Updated for individual patient tracking with age-dependent parameters

library(hesim)
library(data.table)

# Load setup data
load("data/hesim_setup.RData")

cat("Defining model parameters for Phase 2 individual patient model...\n")
cat("Patients:", format(n_patients, big.mark = ","), "| Risk strata:", n_risk_strata, "\n")

# Model timing parameters
n_cycles <- 80      # Lifetime horizon (up to ~50 years from mean age 48.8)
cycle_length <- 1   # 1 year cycles
discount_rate <- 0.03  # 3% annual discount rate

cat("\nModel timing:\n")
cat("- Cycles:", n_cycles, "\n")
cat("- Cycle length:", cycle_length, "year\n")
cat("- Discount rate:", discount_rate * 100, "%\n")

# Intervention parameters (same as Phase 1)
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

# AGE-DEPENDENT PARAMETERS (New for Phase 2)
# ==========================================
# These are critical for individual patient tracking

# Age-dependent mortality rates (non-suicide deaths)
# From US life tables 2017 (reference 23 in paper)
create_age_mortality_table <- function() {
  
  cat("\nCreating age-dependent mortality parameters...\n")
  
  # Age-specific mortality rates per 1000 population (approximate US rates)
  # These will be interpolated for specific ages during simulation
  age_mortality <- data.table(
    age_group = c("18-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85+"),
    age_midpoint = c(21, 29.5, 39.5, 49.5, 59.5, 69.5, 79.5, 90),
    mortality_rate = c(
      0.00087,  # 18-24: 0.87 per 1000
      0.00120,  # 25-34: 1.20 per 1000  
      0.00201,  # 35-44: 2.01 per 1000
      0.00431,  # 45-54: 4.31 per 1000
      0.00965,  # 55-64: 9.65 per 1000
      0.02170,  # 65-74: 21.70 per 1000
      0.05112,  # 75-84: 51.12 per 1000
      0.13137   # 85+: 131.37 per 1000
    )
  )
  
  return(age_mortality)
}

age_mortality_table <- create_age_mortality_table()

# Function to get mortality rate for specific age
get_mortality_rate <- function(age) {
  # Linear interpolation between age groups
  if (age <= 21) return(age_mortality_table$mortality_rate[1])
  if (age >= 90) return(age_mortality_table$mortality_rate[8])
  
  # Find bracketing ages and interpolate
  lower_idx <- max(which(age_mortality_table$age_midpoint <= age))
  upper_idx <- min(which(age_mortality_table$age_midpoint > age))
  
  if (lower_idx == upper_idx) return(age_mortality_table$mortality_rate[lower_idx])
  
  # Linear interpolation
  lower_age <- age_mortality_table$age_midpoint[lower_idx]
  upper_age <- age_mortality_table$age_midpoint[upper_idx]
  lower_rate <- age_mortality_table$mortality_rate[lower_idx]
  upper_rate <- age_mortality_table$mortality_rate[upper_idx]
  
  weight <- (age - lower_age) / (upper_age - lower_age)
  interpolated_rate <- lower_rate + weight * (upper_rate - lower_rate)
  
  return(interpolated_rate)
}

# Test mortality function
cat("Sample mortality rates:\n")
cat("- Age 30:", round(get_mortality_rate(30) * 1000, 2), "per 1000\n")
cat("- Age 50:", round(get_mortality_rate(50) * 1000, 2), "per 1000\n") 
cat("- Age 70:", round(get_mortality_rate(70) * 1000, 2), "per 1000\n")

# Age-dependent cost parameters (2016 USD)
# From Table 1 - background healthcare costs vary by age
cost_params <- list(
  # Suicide attempt costs (age-independent)
  nonfatal_attempt_medical = 10830,    # Medical cost per nonfatal attempt
  nonfatal_attempt_productivity = 17369, # Productivity cost per nonfatal attempt
  fatal_attempt_medical = 4354,        # Medical cost per fatal attempt (age-adjusted)
  fatal_attempt_productivity = 61150,  # Productivity cost per fatal attempt (age-adjusted)
  
  # Background healthcare costs by age group (annual) - KEY FOR INDIVIDUAL TRACKING
  bg_medical_18_44 = 4016,
  bg_medical_45_64 = 7648,
  bg_medical_65plus = 11740,
  
  # Risk assessment cost (evaluation after positive screen)
  evaluation_cost = 76
)

# Function to get background medical cost by age
get_background_cost <- function(age) {
  if (age < 45) return(cost_params$bg_medical_18_44)
  if (age < 65) return(cost_params$bg_medical_45_64)
  return(cost_params$bg_medical_65plus)
}

cat("\nAge-dependent cost parameters (2016 USD):\n")
cat("- Ages 18-44:", cost_params$bg_medical_18_44, "\n")
cat("- Ages 45-64:", cost_params$bg_medical_45_64, "\n")
cat("- Ages 65+:", cost_params$bg_medical_65plus, "\n")
cat("- Nonfatal attempt medical:", cost_params$nonfatal_attempt_medical, "\n")
cat("- Fatal attempt medical:", cost_params$fatal_attempt_medical, "\n")

# Create patient-specific parameter mapping for efficient simulation
create_patient_parameters <- function() {
  
  cat("\nCreating patient-specific parameter mapping...\n")
  cat("Available patient columns:", names(patients), "\n")
  
  # Work with the existing age column name from your setup
  # Your setup script created: patient_id, age, sex, risk_stratum
  age_col <- "age"  # This is what your setup script created
  
  cat("Using age column:", age_col, "\n")
  
  # Create a table linking each patient to their risk stratum parameters
  patient_params <- merge(patients, 
                          risk_strata[, .(stratum_id, baseline_attempt_rate)], 
                          by.x = "risk_stratum", by.y = "stratum_id")
  
  # Add age-dependent parameters using the existing age column
  patient_params[, initial_mortality_rate := sapply(age, get_mortality_rate)]
  patient_params[, initial_bg_cost := sapply(age, get_background_cost)]
  
  # For simulation, we'll need current_age (copy from age initially)
  patient_params[, current_age := age]
  patient_params[, initial_age := age]  # Keep track of starting age too
  
  cat("Patient parameters created for", nrow(patient_params), "patients\n")
  
  return(patient_params)
}

patient_params <- create_patient_parameters()

# Display sample of patient parameters
cat("\nSample patient parameters:\n")
print(head(patient_params[, .(patient_id, age, initial_age, current_age, risk_stratum, baseline_attempt_rate, 
                              initial_mortality_rate, initial_bg_cost)], 10))

# Create transition probability parameters for hesim
# This is more complex for individual patients but follows same structure
create_transition_params <- function() {
  
  cat("\nCreating transition probability parameters for individual patients...\n")
  
  # For hesim, we need transition parameters for each patient-strategy combination
  # This will be handled efficiently in the transition model
  
  # Create base transition data structure
  trans_data <- expand(hesim_dat, by = c("strategies", "patients"))
  
  # Add patient-specific risk parameters
  trans_data <- merge(trans_data[, .(strategy_id, patient_id, strategy_name, age, sex)], 
                      patient_params[, .(patient_id, risk_stratum, baseline_attempt_rate)], 
                      by = "patient_id")
  
  # Add intervention effects
  trans_data[, rr := intervention_params$rr[strategy_name]]
  
  # Calculate adjusted attempt rates
  trans_data[, adjusted_attempt_rate := baseline_attempt_rate * rr]
  
  cat("Transition parameters dataset created with", format(nrow(trans_data), big.mark = ","), "rows\n")
  cat("(", nrow(strategies), "strategies ×", format(n_patients, big.mark = ","), "patients)\n")
  
  return(trans_data)
}

# Create the transition parameters
transition_data <- create_transition_params()

cat("\nSample of transition data:\n")
print(head(transition_data[, .(strategy_name, patient_id, risk_stratum, baseline_attempt_rate, rr, adjusted_attempt_rate)]))

# Save all parameters
cat("\nSaving parameters...\n")
save(
  n_cycles, cycle_length, discount_rate,
  intervention_params, clinical_params, cost_params,
  age_mortality_table, get_mortality_rate, get_background_cost,
  patient_params, transition_data,
  file = "data/hesim_parameters.RData"
)

# cat("\n" + rep("=", 70) + "\n")
cat("PHASE 2 PARAMETERS COMPLETE\n")
# cat(rep("=", 70) + "\n")
cat("Key Features:\n")
cat("✓ Individual patient parameters (", format(n_patients, big.mark = ","), "patients)\n")
cat("✓ Age-dependent mortality rates (increases with age)\n")
cat("✓ Age-stratified background costs (18-44, 45-64, 65+)\n")
cat("✓ Patient-specific risk stratum assignment\n")
cat("✓ Efficient parameter lookup functions\n")

cat("\nParameter Summary:\n")
cat("- Transition combinations:", format(nrow(transition_data), big.mark = ","), "\n")
cat("- Age mortality table:", nrow(age_mortality_table), "age groups\n")
cat("- Cost categories:", length(cost_params), "\n")
cat("- Clinical parameters:", length(clinical_params), "\n")

cat("\nValidation:\n")
cat("- Mean starting age:", round(mean(patient_params$age), 1), "years\n")
cat("- Mean mortality rate:", round(mean(patient_params$initial_mortality_rate) * 1000, 2), "per 1000\n")
cat("- Mean background cost: $", round(mean(patient_params$initial_bg_cost)), "\n")

cat("\nNext: Run 03-transitions.R for age-dependent transition matrices\n")
# cat(rep("=", 70) + "\n")