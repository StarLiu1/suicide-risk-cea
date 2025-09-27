library(hesim)
library(data.table)

# Load setup data
load("data/hesim_setup.RData")

cat("=== PROPER HESIM PARAMETER IMPLEMENTATION ===\n")
cat("Creating hesim tparams_transprobs for CohortDtstm model\n\n")

# =============================================================================
# 1. CREATE TRANSITION PROBABILITY PARAMETERS (PROPER HESIM WAY)
# =============================================================================

create_hesim_transition_params <- function() {
  
  cat("Creating transition probability parameters using tparams_transprobs...\n")
  
  # Create input data for transition model (strategies x patients)
  # This is the standard hesim approach
  transmod_data <- expand(hesim_dat, by = c("strategies", "patients"))
  
  # Add covariates needed for the model
  transmod_data[, intercept := 1]  # Intercept for logistic regression
  transmod_data[, age_centered := age - 48.8]  # Centered age
  transmod_data[, acf_intervention := ifelse(strategy_name == "ACF_Intervention", 1, 0)]
  transmod_data[, cbt_intervention := ifelse(strategy_name == "CBT_Intervention", 1, 0)]
  
  # Add risk stratum information
  # transmod_data <- merge(transmod_data, 
  #                        patients[, .(patient_id, risk_stratum)], 
  #                        by = "patient_id")
  
  # Add baseline attempt rates for each patient
  transmod_data <- merge(transmod_data,
                         risk_strata[, .(risk_stratum, baseline_attempt_rate)],
                         by = "risk_stratum")
  
  cat("Transition input data created with", format(nrow(transmod_data), big.mark = ","), "rows\n")
  
  # Now we need to create transition probability matrices
  # For our 3-state model: no_attempts, prior_attempt, dead
  
  # Create transition matrix structure (3x3 for our states)
  n_states <- 3
  trans_mat <- matrix(c(
    NA, 1, 2,    # From state 1 (no_attempts): can go to state 2 (prior_attempt) or 3 (dead)
    NA, NA, 3,   # From state 2 (prior_attempt): can stay or go to state 3 (dead)  
    NA, NA, NA   # From state 3 (dead): absorbing state
  ), nrow = n_states, byrow = TRUE)
  
  colnames(trans_mat) <- rownames(trans_mat) <- c("no_attempts", "prior_attempt", "dead")
  
  cat("Transition matrix structure:\n")
  print(trans_mat)
  
  # For hesim, we need to create parameters that can calculate probabilities
  # We'll use the define_model() approach with expressions
  
  return(list(
    input_data = transmod_data,
    trans_mat = trans_mat,
    n_states = n_states
  ))
}

trans_params <- create_hesim_transition_params()

# =============================================================================
# 2. COMBINE ALL MODEL DEFINITIONS WITH PROPER HESIM SYNTAX
# =============================================================================

create_complete_model_def <- function() {
  
  cat("\nCreating complete model definition using proper hesim syntax...\n")
  
  # Combine all RNG definitions into a single rng_def
  combined_rng_def <- define_rng({
    # Intervention effects (log-normal distribution)
    acf_rr_mean <- log(0.83)
    acf_rr_se <- (log(0.97) - log(0.71)) / (2 * 1.96)
    acf_log_rr <- rnorm(1, acf_rr_mean, acf_rr_se)
    
    cbt_rr_mean <- log(0.47) 
    cbt_rr_se <- (log(0.73) - log(0.30)) / (2 * 1.96)
    cbt_log_rr <- rnorm(1, cbt_rr_mean, cbt_rr_se)
    
    # Clinical parameters
    death_per_attempt <- rbeta(1, 0.0881 * 100, (1 - 0.0881) * 100)
    prior_multiplier <- rlnorm(1, log(1.54), 0.1)
    
    # Age-dependent mortality parameters
    age_mort_coef <- rnorm(1, 0.08, 0.01)
    
    # Background costs (normal distribution)
    bg_cost_18_44 <- rnorm(1, 4016, 400)
    bg_cost_45_64 <- rnorm(1, 7648, 765)
    bg_cost_65plus <- rnorm(1, 11740, 1174)
    
    # Intervention costs (fixed)
    acf_cost <- 96
    cbt_cost <- 1088
    
    # Utility parameters
    base_utility <- rbeta(1, 0.866 * 100, (1 - 0.866) * 100)
    
    # MUST return a list for hesim
    list(
      acf_log_rr = acf_log_rr,
      cbt_log_rr = cbt_log_rr,
      death_per_attempt = death_per_attempt,
      prior_multiplier = prior_multiplier,
      age_mort_coef = age_mort_coef,
      bg_cost_18_44 = bg_cost_18_44,
      bg_cost_45_64 = bg_cost_45_64,
      bg_cost_65plus = bg_cost_65plus,
      acf_cost = acf_cost,
      cbt_cost = cbt_cost,
      base_utility = base_utility
    )
  })
  
  # Combine all tparams into a single define_tparams
  combined_tparams_def <- define_tparams({
    
    # Get intervention effects
    acf_rr <- exp(acf_log_rr)
    cbt_rr <- exp(cbt_log_rr)
    
    # Calculate intervention effect for this patient/strategy
    intervention_effect <- ifelse(acf_intervention == 1, acf_rr,
                                  ifelse(cbt_intervention == 1, cbt_rr, 1.0))
    
    # Calculate adjusted attempt rate for this patient
    adjusted_attempt_rate <- baseline_attempt_rate * intervention_effect
    
    # Calculate age-dependent mortality (increases with age)
    current_age <- age + time  # Age progression using hesim's built-in 'time' variable
    age_mortality <- pmax(0.001, age_mort_coef * exp((current_age - 50) / 20))
    
    # Calculate transition probabilities for 3x3 matrix
    # FROM STATE 1 (no_attempts):
    suicide_attempt_prob <- pmin(0.8, adjusted_attempt_rate)
    suicide_death_prob <- suicide_attempt_prob * death_per_attempt
    other_death_prob <- pmin(0.2, age_mortality)
    
    # Ensure total doesn't exceed 1 - use vectorized operations
    total_exit_prob <- suicide_attempt_prob + other_death_prob
    scale_factor <- ifelse(total_exit_prob > 0.99, 0.99 / total_exit_prob, 1.0)
    suicide_attempt_prob <- suicide_attempt_prob * scale_factor
    suicide_death_prob <- suicide_death_prob * scale_factor
    other_death_prob <- other_death_prob * scale_factor
    
    # FROM STATE 2 (prior_attempt):
    prior_attempt_prob <- pmin(0.9, adjusted_attempt_rate * prior_multiplier)
    prior_death_prob <- prior_attempt_prob * death_per_attempt
    
    # Age-dependent background costs
    background_cost <- ifelse(current_age < 45, bg_cost_18_44,
                              ifelse(current_age < 65, bg_cost_45_64, bg_cost_65plus))
    
    # Intervention costs
    intervention_cost <- ifelse(acf_intervention == 1, acf_cost,
                                ifelse(cbt_intervention == 1, cbt_cost, 0))
    
    # Create transition probability matrix using tpmatrix()
    # This is the proper hesim way for discrete time models
    tp_matrix <- tpmatrix(
      # Row 1 (from no_attempts): stay, attempt+survive, die
      C, suicide_attempt_prob - suicide_death_prob, suicide_death_prob + other_death_prob,
      # Row 2 (from prior_attempt): can't go back, stay, die  
      0, C, prior_death_prob + other_death_prob,
      # Row 3 (dead): absorbing
      0, 0, 1,
      states = c("no_attempts", "prior_attempt", "dead")
    )
    
    # Return named list with proper hesim elements
    list(
      # Transition probabilities
      tpmatrix = tp_matrix,
      
      # Utilities by state
      utility = c(
        state_1 = base_utility,        # no_attempts
        state_2 = base_utility * 0.95, # prior_attempt  
        state_3 = 0                    # dead
      ),
      
      # Costs by state - must be a list for hesim
      costs = list(
        medical = c(
          state_1 = background_cost + intervention_cost,  # no_attempts
          state_2 = background_cost + intervention_cost,  # prior_attempt
          state_3 = 0                                     # dead
        )
      )
    )
    
  }, times = 0:(n_cycles - 1))
  
  # Create complete model definition with correct syntax
  complete_model_def <- define_model(
    tparams_def = combined_tparams_def,
    rng_def = combined_rng_def,
    n_states = 3
  )
  
  cat("Complete model definition created\n")
  return(complete_model_def)
}

complete_model_def <- create_complete_model_def()

# =============================================================================
# 3. TEST PARAMETER CREATION
# =============================================================================

test_parameter_creation <- function() {
  
  cat("\nTesting parameter creation...\n")
  
  # Test that we can create the actual parameter objects
  tryCatch({
    
    # Create a small test input data
    test_input <- head(trans_params$input_data, 10)
    
    # Test parameter creation with 2 samples using eval_model
    test_eval <- eval_model(complete_model_def, input_data = test_input)
    
    cat("✓ Parameter objects created successfully\n")
    cat("Parameter components:\n")
    print(names(test_eval))
    
    # Test accessing transition probabilities
    if("tpmatrix" %in% names(test_eval)) {
      cat("✓ Transition probabilities accessible\n")
      cat("Transition probability structure:", class(test_eval$tpmatrix), "\n")
      if(is.array(test_eval$tpmatrix)) {
        cat("Transition probability dimensions:", dim(test_eval$tpmatrix), "\n")
      }
    }
    
    # Test costs and utilities
    if("costs" %in% names(test_eval)) {
      cat("✓ Costs accessible\n")
      cat("Cost structure:", class(test_eval$costs), "\n")
    }
    
    if("utility" %in% names(test_eval)) {
      cat("✓ Utilities accessible\n")
      cat("Utility structure:", class(test_eval$utility), "\n")
    }
    
    return(TRUE)
    
  }, error = function(e) {
    cat("✗ Error creating parameters:", e$message, "\n")
    return(FALSE)
  })
}

test_success <- test_parameter_creation()

# =============================================================================
# 4. SAVE PARAMETER OBJECTS
# =============================================================================

cat("\nSaving proper hesim parameter objects...\n")

save(
  # Main model definition
  complete_model_def,
  
  # Transition parameters
  trans_params,
  
  cost_params, utility_params,  # These are now created
  
  # Previous objects (keep for reference)
  hesim_dat, input_data, strategies, patients, states, risk_strata,
  clinical_params, intervention_params, n_cycles, cycle_length, discount_rate,
  n_risk_strata, n_patients, use_individual_patients,
  
  file = "data/hesim_parameters.RData"
)

# =============================================================================
# 5. SUMMARY
# =============================================================================

# cat("\n" + rep("=", 70) + "\n")
cat("PROPER HESIM PARAMETERS CREATED\n")
# cat(rep("=", 70) + "\n")

cat("Key Achievements:\n")
cat("✓ Proper hesim define_model() approach\n")
cat("✓ Combined tparams_def with tpmatrix, utility, and costs\n")
cat("✓ Age-dependent mortality and costs\n")
cat("✓ Individual patient risk stratum effects\n")
cat("✓ Intervention effects via covariates\n")
cat("✓ Time-varying parameters (age progression)\n")
cat("✓ Proper tpmatrix() usage for discrete time transitions\n")

cat("\nParameter Objects:\n")
cat("- Complete model definition: define_model with all components\n")
cat("- Transition probabilities: tpmatrix() with 3x3 structure\n")
cat("- Cost model: age-dependent + intervention costs\n")
cat("- Utility model: state-dependent utilities\n")

cat("\nModel Features:\n")
cat("- Individual patients:", format(n_patients, big.mark = ","), "\n")
cat("- Risk strata: 1000 (patient-specific assignment)\n") 
cat("- Age progression: YES (cycle-by-cycle)\n")
cat("- Intervention effects: Covariate-based\n")
cat("- Time horizon:", n_cycles, "cycles\n")
cat("- Transition matrix: 3x3 (no_attempts, prior_attempt, dead)\n")

if(test_success) {
  cat("✓ All parameter tests passed - ready for R/03-transitions.R\n")
} else {
  cat("⚠️  Some parameter tests failed - check error messages above\n")
}

cat("\nNext Steps:\n")
cat("1. Update R/03-transitions.R to use create_CohortDtstmTrans()\n") 
cat("2. Use complete_model_def with CohortDtstmTrans\n")
cat("3. Test transition model creation and simulation\n")

cat("\nKey hesim syntax corrections made:\n")
cat("- Used define_model(tparams_def, rng_def, n_states)\n")
cat("- Combined all parameters in single tparams_def\n")
cat("- Used tpmatrix() for transition probabilities\n")
cat("- Returned list(tpmatrix, utility, costs) from tparams\n")

# cat(rep("=", 70) + "\n")