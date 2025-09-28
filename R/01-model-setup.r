# 1. Hesim Model Setup - PROPER HESIM IMPLEMENTATION
# File: R/01-model-setup.R
# Full 1000-stratum model using hesim package correctly

library(hesim)
library(data.table)
library(dplyr)
# Load required library
library(readxl)


cat("=== PROPER HESIM IMPLEMENTATION ===\n")
cat("Ross et al. (2021) Suicide Risk Prediction Model\n\n")

# Clear workspace
rm(list = ls())

# Model configuration
n_risk_strata <- 1000  # Full model as per paper
n_strategies <- 1      # Start with one strategy for testing, expand later
use_individual_patients <- TRUE

if (use_individual_patients) {
  # Individual patient model (for age-dependent mortality/costs)
  n_patients <- 100000  # Reasonable size for Phase 2
  cat("Using individual patient model with", n_patients, "patients\n")
} else {
  # Cohort model (simpler, but less accurate for age effects)
  n_patients <- 1
  cat("Using cohort model\n")
}


# =============================================================================
# 7. ADDITIONAL PARAMETERS FOR SIMULATION
# =============================================================================

# Model timing
n_cycles <- 1
cycle_length <- 1  # years
discount_rate <- 0.03

# Clinical parameters from Ross et al.
clinical_params <- list(
  death_per_attempt = 0.0881,  # 8.81% mortality per attempt
  prior_attempt_multiplier = 1.54,  # 54% higher risk for those with prior attempts
  base_utility = 0.866  # EQ-5D utility for primary care population
)

# Intervention parameters
intervention_params <- list(
  rr = c("No_Prediction" = 1.0, "ACF_Intervention" = 0.83, "CBT_Intervention" = 0.47),
  annual_cost = c("No_Prediction" = 0, "ACF_Intervention" = 96, "CBT_Intervention" = 1088),
  uptake = c("No_Prediction" = 1.0, "ACF_Intervention" = 0.994, "CBT_Intervention" = 0.899)
)

# Load and process the life table data
load_mortality_table <- function() {
  
  cat("Loading mortality data from External Data/lifetable.xlsx...\n")
  
  # Read the Excel file
  mortality_data <- read_excel("External Data/lifetable.xlsx")
  
  # Parse age groups and create lookup table
  mortality_lookup <- data.table()
  
  for (i in 1:nrow(mortality_data)) {
    age_group <- mortality_data$`Age Group`[i]
    probability <- mortality_data$Probability[i]
    
    # Parse different age group formats
    if (grepl("–", age_group)) {
      # Format like "99–100"
      ages <- strsplit(age_group, "–")[[1]]
      start_age <- as.numeric(ages[1])
    } else if (grepl("\\+", age_group)) {
      # Format like "100+"
      start_age <- as.numeric(gsub("\\+", "", age_group))
    } else {
      # Assume it's just a number
      start_age <- as.numeric(age_group)
    }
    
    mortality_lookup <- rbind(mortality_lookup, data.table(
      age = start_age,
      mortality_prob = probability
    ))
  }
  
  # Sort by age
  mortality_lookup <- mortality_lookup[order(age)]
  
  cat("Mortality table loaded with", nrow(mortality_lookup), "age groups\n")
  cat("Age range:", min(mortality_lookup$age), "to", max(mortality_lookup$age), "\n")
  
  return(mortality_lookup)
}

# Load the mortality table
mortality_table <- load_mortality_table()

# Age-dependent mortality (for individual patient model)
if (use_individual_patients) {
  # # Simple age-mortality relationship (US life tables approximation)
  # get_age_mortality <- function(age) {
  #   # Approximate annual mortality rates by age
  #   ifelse(age < 25, 0.001,
  #          ifelse(age < 35, 0.0015,
  #                 ifelse(age < 45, 0.002,
  #                        ifelse(age < 55, 0.005,
  #                               ifelse(age < 65, 0.01,
  #                                      ifelse(age < 75, 0.025,
  #                                             ifelse(age < 85, 0.06, 0.15)))))))
  # }
  # Replace get_age_mortality function entirely
  get_age_mortality <- function(age) {
    
    # Floor age to integer
    age_int <- floor(age)
    
    # Handle edge cases
    if (age_int < 0) return(0)
    if (age_int >= max(mortality_table$age)) {
      # Use highest available age group for very old ages
      return(mortality_table[age == max(mortality_table$age)]$mortality_prob)
    }
    
    # Find the appropriate age group
    # Use the highest age that is <= input age
    applicable_ages <- mortality_table[age <= age_int]
    
    if (nrow(applicable_ages) == 0) {
      # Age is below minimum in table, use first entry
      return(mortality_table[1]$mortality_prob)
    } else {
      # Use the most recent applicable age group
      return(applicable_ages[age == max(applicable_ages$age)]$mortality_prob)
    }
  }
  get_age_mortality <- Vectorize(get_age_mortality)
  
}

suicide_mortality_data <- read_excel("External Data/age_specific_suicide.xlsx")

# Create lookup table
suicide_mortality_table <- data.table(
  age = suicide_mortality_data$Age,
  suicide_rate = suicide_mortality_data$`% of Total Deaths`  # Already numeric
)

get_suicide_percentage <- function(age) {
  age_int <- floor(age)
  
  # Handle edge cases
  if (age_int < min(suicide_mortality_table$age)) {
    return(suicide_mortality_table[1]$suicide_rate)
  }
  if (age_int > max(suicide_mortality_table$age)) {
    return(suicide_mortality_table[nrow(suicide_mortality_table)]$suicide_rate)
  }
  
  # Find closest age
  closest_age <- suicide_mortality_table[age <= age_int][which.max(age)]
  return(closest_age$suicide_rate)
}
get_suicide_percentage <- Vectorize(get_suicide_percentage)

get_suicide_mortality_probability <- function(age) {
  all_cause_prob <- get_age_mortality(age)  # e.g., 0.005 probability
  suicide_percentage <- get_suicide_percentage(age)  # e.g., 0.02 (2% of deaths)
  
  # Among those who die, suicide_percentage die from suicide
  # So: P(suicide death) = P(any death) × P(suicide | death)
  suicide_prob <- all_cause_prob * suicide_percentage
  
  return(suicide_prob)
}

get_background_mortality <- function(age) {
  all_cause_prob <- get_age_mortality(age)  # Your existing function
  suicide_prob <- get_suicide_mortality_probability(age)
  
  # Subtract suicide deaths to get non-suicide background mortality
  background_prob <- all_cause_prob - suicide_prob
  
  # Ensure non-negative
  return(pmax(0, background_prob))
}
get_background_mortality <- Vectorize(get_background_mortality)

# =============================================================================
# 1. STRATEGIES (INTERVENTIONS)
# =============================================================================

strategies <- data.table(
  strategy_id = 1:3,
  strategy_name = c("No_Prediction", "ACF_Intervention", "CBT_Intervention")
)

cat("\nStrategies defined:\n")
print(strategies)

# =============================================================================
# 2. PATIENTS 
# =============================================================================

if (use_individual_patients) {
  # Generate individual patients with age variation and risk stratum assignment
  set.seed(12345)  # For reproducibility
  
  patients <- data.table(
    patient_id = 1:n_patients,
    # Age distribution from paper (mean 48.8, SD 17.2)
    # age = rnorm(n_patients, 48.8, 17.2),
    age = pmax(0, pmin(100, rnorm(n_patients, 48.8, 17.2))),
    # Risk stratum (1-1000, uniform distribution)
    risk_stratum = sample(1:n_risk_strata, n_patients, replace = TRUE)
  )
  
  cat("\nPatients created:\n")
  cat("- Count:", nrow(patients), "\n")
  cat("- Age: mean =", round(mean(patients$age), 2), ", range =", 
      round(min(patients$age)), "-", round(max(patients$age)), "\n")
  
} else {
  # Single representative patient for cohort model
  patients <- data.table(
    patient_id = 1,
    age = 48.8,  # Mean age
    risk_stratum = NA  # Will be handled differently in cohort model
  )
}

# =============================================================================
# 3. HEALTH STATES (HESIM FORMAT)
# =============================================================================

# For hesim, we need a simpler state structure
# Each risk stratum has 3 states: no_attempts, prior_attempt, dead

if (use_individual_patients) {
  # With individual patients, states are simpler since risk stratum is patient attribute
  states <- data.table(
    state_id = 1:3,
    state_name = c("no_attempts", "prior_attempt", "dead")
  )
} else {
  # For cohort model, we need separate states for each risk stratum
  states <- data.table()
  for (stratum in 1:n_risk_strata) {
    stratum_states <- data.table(
      state_id = ((stratum-1)*3 + 1):((stratum-1)*3 + 3),
      state_name = paste0("s", stratum, "_", c("no_attempts", "prior_attempt", "dead")),
      risk_stratum = stratum
    )
    states <- rbind(states, stratum_states)
  }
}

cat("\nHealth states defined:\n")
cat("- Total states:", nrow(states), "\n")
print(head(states))

# =============================================================================
# 4. RISK STRATA PARAMETERS (FROM ROSS ET AL.)
# =============================================================================

create_ross_risk_distribution <- function() {
  cat("\nCreating Ross et al. risk distribution (1000 strata)...\n")
  
  # Exact rates from Table 1 (converted from percentages to decimals)
  percentile_rates <- data.table(
    percentile_range = c("0-90", "90-95", "95-99", ">99"),
    rate_percent = c(0.008, 0.213, 1.351, 15.437),
    rate_decimal = c(0.008, 0.213, 1.351, 15.437) / 100
  )
  
  print(percentile_rates)
  
  # Create 1000 strata
  risk_strata <- data.table(
    risk_stratum = 1:n_risk_strata,
    percentile = (1:n_risk_strata - 0.5) / n_risk_strata  # Midpoint percentiles
  )
  
  # Assign rates based on percentile brackets
  risk_strata[, baseline_attempt_rate := case_when(
    percentile <= 0.90 ~ 0.00008,   # 0-90th percentile: 0.008% 
    percentile <= 0.95 ~ 0.00213,   # 90-95th percentile: 0.213%
    percentile <= 0.99 ~ 0.01351,   # 95-99th percentile: 1.351%
    TRUE ~ 0.15437                  # >99th percentile: 15.437%
  )]
  
  # Population weights (each stratum is 1/1000)
  risk_strata[, population_weight := 1/n_risk_strata]
  
  # Validate against target rates
  overall_attempt_rate <- sum(risk_strata$baseline_attempt_rate * risk_strata$population_weight)
  target_rate <- 175 / 100000  # 175 per 100,000 person-years
  
  cat("Expected attempt rate:", round(overall_attempt_rate * 100000, 1), "per 100,000\n")
  cat("Target attempt rate: 175 per 100,000\n")
  
  # Check death rate (target: 15 per 100k person-years)  
  expected_death_rate <- overall_attempt_rate * clinical_params$death_per_attempt
  cat("Expected death rate:", round(expected_death_rate * 100000, 1), "per 100k (target: 15)\n")
  target_death_rate <- 15 / 100000  # 175 per 100,000 person-years
  
  # Apply calibration factor if needed
  calibration_factor <- target_rate / overall_attempt_rate
  risk_strata[, baseline_attempt_rate := baseline_attempt_rate * calibration_factor]
  
  final_rate <- sum(risk_strata$baseline_attempt_rate * risk_strata$population_weight)
  cat("Calibrated attempt rate:", round(final_rate * 100000, 1), "per 100,000\n")
  
  expected_death_rate <- final_rate * clinical_params$death_per_attempt
  
  cat("Calibrated death rate:", round(expected_death_rate * 100000, 1), "per 100k (target: 15)\n")

  return(risk_strata)
}

risk_strata <- create_ross_risk_distribution()

# =============================================================================
# 5. CREATE HESIM DATA OBJECT
# =============================================================================

cat("\nCreating hesim_data object...\n")

# This is the core hesim data structure
hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients,
  states = states
)

cat("Hesim data object created successfully:\n")
cat("- Strategies:", nrow(hesim_dat$strategies), "\n")
cat("- Patients:", nrow(hesim_dat$patients), "\n")
cat("- States:", nrow(hesim_dat$states), "\n")

# =============================================================================
# 6. CREATE EXPANDED INPUT DATA FOR HESIM MODELS
# =============================================================================

# This creates all strategy × patient × state combinations for modeling
input_data <- expand(hesim_dat)
cat("Input data expanded to", nrow(input_data), "rows\n")

# Add risk stratum information to input data
if (use_individual_patients) {
  # Risk stratum comes from patient data
  input_data <- merge(input_data[, .(patient_id, strategy_id, strategy_name, age)], 
                      patients[, .(patient_id, risk_stratum)], 
                      by = "patient_id")
} else {
  # Risk stratum comes from state data
  input_data <- merge(input_data[, .(patient_id, strategy_id, strategy_name, age)], 
                      states[, .(state_id, risk_stratum)], 
                      by = "state_id")
}

# Add baseline attempt rates
input_data <- merge(input_data[, .(patient_id, strategy_id, strategy_name, age, risk_stratum)], 
                    risk_strata[, .(risk_stratum, baseline_attempt_rate)],
                    by = "risk_stratum")

cat("Input data enriched with risk parameters\n")


# =============================================================================
# 8. SAVE SETUP DATA
# =============================================================================

cat("\nSaving setup data...\n")

if (!dir.exists("data")) dir.create("data")

save(
  # Core hesim objects
  hesim_dat, input_data,
  
  # Model structure
  strategies, patients, states, risk_strata,
  
  # Parameters
  clinical_params, intervention_params,
  n_cycles, cycle_length, discount_rate,
  
  # Configuration
  n_risk_strata, n_patients, use_individual_patients,
  
  # Functions (if individual patient model)
  list = if(use_individual_patients) c(ls(), "get_age_mortality") else ls(),
  
  file = "data/hesim_setup.RData"
)

# cat("\n" + strrep("=", 70) + "\n")
cat("HESIM SETUP COMPLETE\n")
# cat(strrep("=", 70) + "\n")

cat("Configuration:\n")
cat("- Model type:", if(use_individual_patients) "Individual Patient" else "Cohort", "\n")
cat("- Risk strata:", n_risk_strata, "\n")
cat("- Patients:", n_patients, "\n")
cat("- States:", nrow(states), "\n")
cat("- Strategies:", nrow(strategies), "\n")
cat("- Time horizon:", n_cycles, "cycles\n")

cat("\nKey Features:\n")
cat("✓ Proper hesim data structure\n")
cat("✓ Ross et al. risk distribution (1000 strata)\n")
cat("✓ Individual patient tracking with age\n") 
cat("✓ All parameters from paper\n")
cat("✓ Ready for hesim transition/cost/utility models\n")

cat("\nNext: Run R/02-parameters.R to set up hesim model parameters properly\n")
# cat(strrep("=", 70) + "\n")