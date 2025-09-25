# Corrected Hesim Simulation - Individual Patient Model
# File: R/corrected-05-simulation.R

library(hesim)
library(data.table)

# Load all components
load("data/hesim_setup.RData")
load("data/hesim_parameters.RData")
load("data/corrected_hesim_transitions.RData")
load("data/hesim_costs_utilities.RData")

cat("Running CORRECTED individual patient simulation...\n")
cat("Model: Individual patient tracking with age progression\n\n")

# SIMULATION PARAMETERS
# =====================

n_cycles_sim <- 20  # Start with shorter horizon for testing
n_samples_sim <- 1  # Deterministic run

cat("Simulation setup:\n")
cat("- Patients:", format(n_patients, big.mark = ","), "\n")
cat("- Strategies:", nrow(strategies), "\n")
cat("- Cycles:", n_cycles_sim, "\n")
cat("- Samples:", n_samples_sim, "\n\n")

# RUN STATE PROBABILITY SIMULATION
# =================================

cat("=== RUNNING STATE PROBABILITY SIMULATION ===\n")

# Use our corrected transition model
transition_model$sim_stateprobs(n_cycles = n_cycles_sim)

# Get results
state_probs <- transition_model$stateprobs_

cat("State probabilities simulated successfully!\n")
cat("Results dimensions:", nrow(state_probs), "rows\n")

# ANALYZE RESULTS
# ===============

analyze_simulation_results <- function(state_probs) {
  
  cat("\n=== ANALYZING SIMULATION RESULTS ===\n")
  
  # Calculate key outcomes by strategy
  results_by_strategy <- state_probs[, {
    
    # Calculate population in each state at each time point
    pop_no_attempts <- sum(prob[state_id == 1])
    pop_prior_attempt <- sum(prob[state_id == 2]) 
    pop_dead <- sum(prob[state_id == 3])
    
    list(
      pop_no_attempts = pop_no_attempts,
      pop_prior_attempt = pop_prior_attempt,
      pop_dead = pop_dead,
      total_pop = pop_no_attempts + pop_prior_attempt + pop_dead
    )
    
  }, by = .(sample, strategy_id, t)]
  
  # Add strategy names
  results_by_strategy <- merge(results_by_strategy,
                               strategies[, .(strategy_id, strategy_name)],
                               by = "strategy_id")
  
  # Calculate cumulative outcomes over time
  final_results <- results_by_strategy[t == max(t), {
    
    # Calculate rates per 100,000 person-years
    person_years <- n_patients * n_cycles_sim
    
    # Attempts = people who transitioned to prior_attempt state
    total_attempts <- pop_prior_attempt + pop_dead  # Approximation
    attempt_rate_per_100k <- (total_attempts / person_years) * 100000
    
    # Deaths = people in dead state
    death_rate_per_100k <- (pop_dead / person_years) * 100000
    
    list(
      final_alive = pop_no_attempts + pop_prior_attempt,
      final_dead = pop_dead,
      total_attempts = total_attempts,
      attempt_rate_per_100k = attempt_rate_per_100k,
      death_rate_per_100k = death_rate_per_100k
    )
    
  }, by = .(strategy_name)]
  
  return(list(
    by_strategy = results_by_strategy,
    final_results = final_results
  ))
}

# Analyze results
analysis_results <- analyze_simulation_results(state_probs)

# Display final results
cat("\n=== FINAL RESULTS BY STRATEGY ===\n")
print(analysis_results$final_results)

# CALCULATE COST AND UTILITY OUTCOMES
# ====================================

calculate_economic_outcomes <- function(state_probs) {
  
  cat("\n=== CALCULATING ECONOMIC OUTCOMES ===\n")
  
  # For each patient-strategy-cycle combination, calculate costs and utilities
  economic_results <- state_probs[, {
    
    # Get strategy information
    strat_name <- strategies$strategy_name[strategy_id]
    
    # Calculate costs for this patient in this cycle
    # (Simplified - in full model would use patient-specific age-dependent costs)
    
    # Base costs
    intervention_cost <- intervention_params$annual_cost[[strat_name]]
    background_cost <- 5000  # Approximate average background cost
    
    # State-specific costs and utilities
    if (state_id == 1) {  # No attempts
      state_cost <- background_cost + intervention_cost
      state_utility <- clinical_params$base_utility
    } else if (state_id == 2) {  # Prior attempt
      state_cost <- background_cost + intervention_cost + cost_params$nonfatal_attempt_medical
      state_utility <- clinical_params$base_utility * 0.95
    } else {  # Dead
      state_cost <- 0
      state_utility <- 0
    }
    
    # Apply discounting
    discount_factor <- 1 / (1 + discount_rate)^t
    discounted_cost <- state_cost * discount_factor * prob
    discounted_utility <- state_utility * discount_factor * prob
    
    list(
      discounted_cost = discounted_cost,
      discounted_utility = discounted_utility
    )
    
  }, by = .(sample, strategy_id, patient_id, state_id, t)]
  
  # Sum across patients and time for each strategy
  strategy_totals <- economic_results[, {
    list(
      total_cost = sum(discounted_cost),
      total_utility = sum(discounted_utility)
    )
  }, by = .(sample, strategy_id)]
  
  # Add strategy names and calculate per-patient values
  strategy_totals[, strategy_name := strategies$strategy_name[strategy_id]]
  strategy_totals[, cost_per_patient := total_cost / n_patients]
  strategy_totals[, utility_per_patient := total_utility / n_patients]
  
  return(strategy_totals)
}

# Calculate economic outcomes
economic_outcomes <- calculate_economic_outcomes(state_probs)

cat("=== ECONOMIC OUTCOMES BY STRATEGY ===\n")
print(economic_outcomes[, .(strategy_name, cost_per_patient, utility_per_patient)])

# VALIDATION AGAINST ROSS ET AL. TARGETS
# =======================================

validate_against_targets <- function(final_results) {
  
  cat("\n=== VALIDATION AGAINST ROSS ET AL. TARGETS ===\n")
  
  # Get baseline results (No_Prediction strategy)
  baseline <- final_results[strategy_name == "No_Prediction"]
  
  if (nrow(baseline) == 0) {
    cat("Warning: No baseline strategy found\n")
    return(NULL)
  }
  
  cat("Ross et al. targets:\n")
  cat("- Suicide attempts: 175 per 100,000 person-years\n") 
  cat("- Suicide deaths: 15 per 100,000 person-years\n\n")
  
  cat("Our simulation results (baseline):\n")
  cat(sprintf("- Suicide attempts: %.1f per 100,000 person-years\n", 
              baseline$attempt_rate_per_100k))
  cat(sprintf("- Suicide deaths: %.1f per 100,000 person-years\n", 
              baseline$death_rate_per_100k))
  
  # Calculate validation ratios
  attempt_ratio <- baseline$attempt_rate_per_100k / 175
  death_ratio <- baseline$death_rate_per_100k / 15
  
  cat(sprintf("\nValidation ratios:\n"))
  cat(sprintf("- Attempt rate ratio: %.3f (target: 1.000)\n", attempt_ratio))
  cat(sprintf("- Death rate ratio: %.3f (target: 1.000)\n", death_ratio))
  
  # Assessment
  if (attempt_ratio >= 0.5 && attempt_ratio <= 2.0) {
    cat("✓ Attempt rates within reasonable range\n")
  } else {
    cat("⚠️  Attempt rates outside expected range\n")
  }
  
  if (death_ratio >= 0.5 && death_ratio <= 2.0) {
    cat("✓ Death rates within reasonable range\n")
  } else {
    cat("⚠️  Death rates outside expected range\n")
  }
  
  return(data.table(
    attempt_ratio = attempt_ratio,
    death_ratio = death_ratio,
    validation_status = ifelse(
      attempt_ratio >= 0.5 && attempt_ratio <= 2.0 && death_ratio >= 0.5 && death_ratio <= 2.0,
      "PASS", "NEEDS_ADJUSTMENT"
    )
  ))
}

# Validate results
validation_results <- validate_against_targets(analysis_results$final_results)

# SAVE RESULTS
# ============

cat("\nSaving corrected simulation results...\n")

if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/results")) dir.create("output/results")

save(
  state_probs,
  analysis_results,
  economic_outcomes,
  validation_results,
  file = "output/results/corrected_individual_patient_results.RData"
)

# CREATE SUMMARY REPORT
# ======================

create_summary_report <- function() {
  
  # cat("\n" + rep("=", 80) + "\n")
  cat("CORRECTED INDIVIDUAL PATIENT SIMULATION SUMMARY\n")
  # cat(rep("=", 80) + "\n")
  
  cat("Model Configuration:\n")
  cat("- Individual patients:", format(n_patients, big.mark = ","), "\n")
  cat("- Risk strata: 1000 (individual assignment)\n")
  cat("- Age progression: YES (cycle-by-cycle)\n")
  cat("- Transition matrices: Patient-specific\n")
  cat("- Time horizon:", n_cycles_sim, "cycles\n\n")
  
  cat("Key Results:\n")
  for (i in 1:nrow(analysis_results$final_results)) {
    result <- analysis_results$final_results[i]
    cat(sprintf("- %s: %.1f attempts/100k, %.1f deaths/100k\n",
                result$strategy_name,
                result$attempt_rate_per_100k,
                result$death_rate_per_100k))
  }
  
  cat("\nValidation Status:\n")
  if (!is.null(validation_results)) {
    cat("- Overall:", validation_results$validation_status, "\n")
    cat("- Attempt rate calibration:", ifelse(validation_results$attempt_ratio > 0.8 & validation_results$attempt_ratio < 1.2, "GOOD", "NEEDS_WORK"), "\n")
    cat("- Death rate calibration:", ifelse(validation_results$death_ratio > 0.8 & validation_results$death_ratio < 1.2, "GOOD", "NEEDS_WORK"), "\n")
  }
  
  cat("\nNext Steps:\n")
  if (!is.null(validation_results) && validation_results$validation_status == "PASS") {
    cat("✓ Model validation successful - ready to extend to full time horizon\n")
    cat("✓ Can proceed with probabilistic sensitivity analysis\n")
    cat("✓ Ready for cost-effectiveness analysis\n")
  } else {
    cat("⚠️  Model calibration needs adjustment\n")
    cat("- Check risk distribution parameters\n")
    cat("- Verify transition matrix calculations\n")
    cat("- Consider scaling factors for target rates\n")
  }
  
  # cat(rep("=", 80) + "\n")
}

create_summary_report()

cat("\nCorrected individual patient simulation completed!\n")
cat("Results saved to: output/results/corrected_individual_patient_results.RData\n")