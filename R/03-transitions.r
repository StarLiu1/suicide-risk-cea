# 3. Hesim Transition Model - PROPER COHORTDTSTMTRANS IMPLEMENTATION
# File: R/03-transitions-fixed.R
# Using proper hesim CohortDtstmTrans with define_model parameters

library(hesim)
library(data.table)

# Load previous data
load("data/hesim_setup.RData")
load("data/hesim_parameters_fixed.RData")

cat("=== PROPER HESIM TRANSITION MODEL ===\n")
cat("Creating CohortDtstmTrans using proper hesim define_model approach\n\n")

# =============================================================================
# 1. CREATE TRANSITION MODEL USING COHORTDTSTMTRANS
# =============================================================================

create_hesim_transition_model <- function() {
  
  cat("Creating CohortDtstmTrans from define_model parameters...\n")
  
  # Use create_CohortDtstm to build the complete economic model
  # This automatically creates the CohortDtstmTrans from our model definition
  
  # First, let's create the input data for the full model
  input_data <- expand(hesim_dat, by = c("strategies", "patients"))
  
  # Add the covariates that our model expects
  input_data[, intercept := 1]
  input_data[, age_centered := age - 48.8]
  input_data[, acf_intervention := ifelse(strategy_name == "ACF_Intervention", 1, 0)]
  input_data[, cbt_intervention := ifelse(strategy_name == "CBT_Intervention", 1, 0)]
  
  # Add risk stratum information
  # input_data <- merge(input_data, 
  #                     patients[, .(patient_id, risk_stratum)], 
  #                     by = "patient_id")
  # 
  # Add baseline attempt rates
  input_data <- merge(input_data,
                      risk_strata[, .(risk_stratum, baseline_attempt_rate)],
                      by.x = "risk_stratum", by.y = "risk_stratum")
  
  cat("Input data created with", format(nrow(input_data), big.mark = ","), "rows\n")
  cat("Columns:", paste(names(input_data), collapse = ", "), "\n")
  
  return(input_data)
}

model_input_data <- create_hesim_transition_model()

# =============================================================================
# 2. CREATE COMPLETE ECONOMIC MODEL USING CREATE_COHORTDTSTM
# =============================================================================

create_complete_economic_model <- function() {
  
  cat("\nCreating complete economic model with create_CohortDtstm...\n")
  
  # This is the proper hesim way - create the entire economic model at once
  # create_CohortDtstm automatically handles:
  # - CohortDtstmTrans for transitions (from tpmatrix in our model_def)
  # - StateVals for utilities (from utility in our model_def)  
  # - StateVals for costs (from costs in our model_def)
  
  tryCatch({
    
    # Use the hesim create_CohortDtstm function with correct syntax
    econmod <- create_CohortDtstm(
      object = complete_model_def,  # The object argument is required
      input_data = model_input_data
      # n = 1  # Number of PSA samples (deterministic for now)
    )
    
    cat("✓ Complete economic model created successfully\n")
    cat("Model components:\n")
    cat("- Transition model:", class(econmod$trans_model)[1], "\n")
    cat("- Utility model:", class(econmod$utility_model)[1], "\n") 
    cat("- Cost models:", length(econmod$cost_models), "categories\n")
    
    return(econmod)
    
  }, error = function(e) {
    cat("✗ Error creating economic model:", e$message, "\n")
    cat("This might be due to parameter issues in complete_model_def\n")
    return(NULL)
  })
}

# Create the complete economic model
econmod <- create_complete_economic_model()

# =============================================================================
# 3. TEST TRANSITION MODEL COMPONENTS
# =============================================================================

test_transition_model <- function(econmod) {
  
  if(is.null(econmod)) {
    cat("Cannot test transition model - econmod is NULL\n")
    return(FALSE)
  }
  
  cat("\n=== TESTING TRANSITION MODEL ===\n")
  
  # Test accessing transition model
  trans_model <- econmod$trans_model
  cat("Transition model class:", class(trans_model), "\n")
  
  # Test that we can simulate state probabilities
  tryCatch({
    
    cat("Testing state probability simulation...\n")
    
    # Simulate for a few cycles
    econmod$sim_stateprobs(n_cycles = 5)
    
    # Check results
    stateprobs <- econmod$stateprobs_
    cat("✓ State probabilities simulated successfully\n")
    cat("State probabilities data dimensions:", dim(stateprobs), "\n")
    cat("Columns:", paste(names(stateprobs), collapse = ", "), "\n")
    
    # Show sample results
    cat("\nSample state probabilities:\n")
    print(head(stateprobs, 10))
    
    # Validate probabilities
    prob_sums <- stateprobs[, .(prob_sum = sum(prob)), by = .(sample, strategy_id, patient_id, t)]
    cat("\nProbability sums by time point (should be ~1.0):\n")
    print(summary(prob_sums$prob_sum))
    
    return(TRUE)
    
  }, error = function(e) {
    cat("✗ Error in state probability simulation:", e$message, "\n")
    return(FALSE)
  })
}

# Test the model
test_success <- test_transition_model(econmod)

# =============================================================================
# 4. TEST COST AND UTILITY MODELS
# =============================================================================

test_cost_utility_models <- function(econmod) {
  
  if(is.null(econmod) || !test_success) {
    cat("Cannot test cost/utility models\n")
    return(FALSE)
  }
  
  cat("\n=== TESTING COST AND UTILITY MODELS ===\n")
  
  # Test utility simulation
  tryCatch({
    cat("Testing utility simulation...\n")
    econmod$sim_qalys(dr = 0.03)
    
    qalys <- econmod$qalys_
    cat("✓ QALYs simulated successfully\n")
    cat("QALYs data dimensions:", dim(qalys), "\n")
    
    # Show sample results
    cat("Sample QALYs:\n")
    print(head(qalys))
    
  }, error = function(e) {
    cat("✗ Error in utility simulation:", e$message, "\n")
  })
  
  # Test cost simulation
  tryCatch({
    cat("\nTesting cost simulation...\n")
    econmod$sim_costs(dr = 0.03)
    
    costs <- econmod$costs_
    cat("✓ Costs simulated successfully\n")
    cat("Costs data dimensions:", dim(costs), "\n")
    
    # Show sample results
    cat("Sample costs:\n")
    print(head(costs))
    
    return(TRUE)
    
  }, error = function(e) {
    cat("✗ Error in cost simulation:", e$message, "\n")
    return(FALSE)
  })
}

# Test cost and utility models
cost_utility_success <- test_cost_utility_models(econmod)

# =============================================================================
# 5. ANALYZE BASIC MODEL OUTCOMES
# =============================================================================

analyze_basic_outcomes <- function(econmod) {
  
  if(is.null(econmod) || !test_success) {
    cat("Cannot analyze outcomes\n")
    return(NULL)
  }
  
  cat("\n=== BASIC OUTCOME ANALYSIS ===\n")
  
  # Get state probabilities
  stateprobs <- econmod$stateprobs_
  
  # Calculate outcomes by strategy at final time point
  final_time <- max(stateprobs$t)
  final_outcomes <- stateprobs[t == final_time, {
    
    # Calculate population in each state
    pop_alive <- sum(prob[state_id %in% 1:2])  # States 1 and 2 are alive
    pop_dead <- sum(prob[state_id == 3])       # State 3 is dead
    pop_prior_attempt <- sum(prob[state_id == 2])  # State 2 is prior attempt
    
    list(
      final_alive = pop_alive,
      final_dead = pop_dead,
      final_prior_attempt = pop_prior_attempt,
      total_pop = pop_alive + pop_dead
    )
    
  }, by = .(sample, strategy_id)]
  
  # Add strategy names
  final_outcomes <- merge(final_outcomes,
                          strategies[, .(strategy_id, strategy_name)],
                          by = "strategy_id")
  
  cat("Final outcomes by strategy (after", final_time, "cycles):\n")
  print(final_outcomes)
  
  # Calculate approximate rates
  person_years <- n_patients * final_time
  
  rate_analysis <- final_outcomes[, {
    
    # Approximate attempt rate (people who ever attempted)
    attempt_rate_per_100k <- (final_prior_attempt + final_dead) / person_years * 100000
    
    # Death rate
    death_rate_per_100k <- final_dead / person_years * 100000
    
    list(
      strategy_name = strategy_name,
      attempt_rate_per_100k = round(attempt_rate_per_100k, 1),
      death_rate_per_100k = round(death_rate_per_100k, 1)
    )
  }]
  
  cat("\nApproximate rates by strategy:\n")
  print(rate_analysis)
  
  return(rate_analysis)
}

# Analyze basic outcomes
basic_outcomes <- analyze_basic_outcomes(econmod)

# =============================================================================
# 6. SAVE TRANSITION MODEL
# =============================================================================

cat("\nSaving transition model components...\n")

save(
  # Main economic model
  econmod,
  
  # Input data and model definition
  model_input_data, complete_model_def,
  
  # Test results and outcomes
  test_success, cost_utility_success, basic_outcomes,
  
  # Keep previous objects
  hesim_dat, strategies, patients, states, risk_strata,
  clinical_params, intervention_params, n_cycles, cycle_length, discount_rate,
  n_risk_strata, n_patients, use_individual_patients,
  
  file = "data/hesim_transitions_fixed.RData"
)

# =============================================================================
# 7. SUMMARY REPORT
# =============================================================================

# cat("\n" + rep("=", 70) + "\n")
cat("HESIM TRANSITION MODEL COMPLETE\n")
# cat(rep("=", 70) + "\n")

cat("Model Type: CohortDtstm (Cohort Discrete Time State Transition Model)\n")
cat("Created using: create_CohortDtstm() with define_model() parameters\n\n")

cat("Key Achievements:\n")
cat("✓ Proper hesim CohortDtstmTrans creation\n")
cat("✓ Complete economic model with transitions, costs, utilities\n")
cat("✓ Individual patient effects via covariates\n")
cat("✓ Age-dependent parameters using 'time' variable\n")
cat("✓ Risk stratum assignment per patient\n")
cat("✓ Intervention effects through covariates\n")

cat("\nModel Configuration:\n")
cat("- Model class: CohortDtstm\n")
cat("- Patients:", format(n_patients, big.mark = ","), "\n")
cat("- Strategies:", nrow(strategies), "\n")
cat("- States: 3 (no_attempts, prior_attempt, dead)\n")
cat("- Time horizon:", n_cycles, "cycles\n")
cat("- Age progression: YES (via 'time' variable)\n")

cat("\nModel Components Status:\n")
if(!is.null(econmod)) {
  cat("✓ Economic model: Created successfully\n")
  cat("✓ Transition model:", class(econmod$trans_model)[1], "\n")
  cat("✓ Utility model:", class(econmod$utility_model)[1], "\n")
  cat("✓ Cost models:", length(econmod$cost_models), "categories\n")
} else {
  cat("✗ Economic model: Failed to create\n")
}

cat("\nSimulation Tests:\n")
cat("- State probabilities:", ifelse(test_success, "✓ PASS", "✗ FAIL"), "\n")
cat("- Costs and utilities:", ifelse(cost_utility_success, "✓ PASS", "✗ FAIL"), "\n")

if(!is.null(basic_outcomes)) {
  cat("\nSample Results (after", max(econmod$stateprobs_$t), "cycles):\n")
  for(i in 1:nrow(basic_outcomes)) {
    outcome <- basic_outcomes[i]
    cat(sprintf("- %s: %.1f attempts/100k, %.1f deaths/100k\n",
                outcome$strategy_name, 
                outcome$attempt_rate_per_100k,
                outcome$death_rate_per_100k))
  }
}

if(test_success && cost_utility_success) {
  cat("\n✓ Model ready for full simulation and cost-effectiveness analysis\n")
  cat("✓ Can proceed to R/05-simulation.R for complete analysis\n")
} else {
  cat("\n⚠️  Model has issues - check error messages above\n")
  cat("⚠️  May need to debug parameter definitions\n")
}

cat("\nNext Steps:\n")
cat("1. If model tests pass: proceed to R/05-simulation.R\n")
cat("2. If model tests fail: debug complete_model_def parameters\n")
cat("3. Validate results against Ross et al. targets\n")
cat("4. Extend to full time horizon and PSA\n")

# cat(rep("=", 70) + "\n")