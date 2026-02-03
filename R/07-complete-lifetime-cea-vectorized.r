# OPTIMIZED Vectorized Complete Lifetime CEA Simulation
# File: R/07-complete-lifetime-cea-vectorized-optimized.R
# With progress reporting and performance optimization

library(hesim)
library(data.table)
library(ggplot2)

load("data/hesim_costs_utilities.RData") #top 5%
# load("data/mat_cal/hesim_costs_utilities.RData") # top 1%

cat("=== OPTIMIZED VECTORIZED LIFETIME CEA SIMULATION ===\n")
cat("Ross et al. (2021) Suicide Risk Prediction Model\n\n")

# =============================================================================
# 1. SIMULATION CONFIGURATION
# =============================================================================

n_cycles_lifetime <- 60
n_samples_psa <- 1

cat("Lifetime simulation configuration:\n")
cat("- Patients:", format(n_patients, big.mark = ","), "\n")
cat("- Risk strata:", n_risk_strata, "\n")
cat("- Time horizon:", n_cycles_lifetime, "cycles\n\n")

# =============================================================================
# 2. OPTIMIZED VECTORIZED SIMULATION ENGINE
# =============================================================================

run_lifetime_simulation_optimized <- function(strategy_name, patients_with_pred, verbose = TRUE) {
  
  if (verbose) cat("\nSimulating strategy:", strategy_name, "...\n")
  
  start_time <- Sys.time()
  
  # Get intervention parameters
  intervention_rr <- intervention_params$rr[[strategy_name]]
  
  # =============================================================================
  # STEP 1: CREATE PATIENT-LEVEL DATA (small, fast)
  # =============================================================================
  
  if (verbose) cat("  [1/5] Preparing patient data... ")
  step_start <- Sys.time()
  
  # Create enhanced patient table with all needed attributes
  patient_data <- copy(patients_with_pred)
  
  # Add risk parameters (merge once on small table)
  patient_data <- merge(
    patient_data,
    risk_strata[, .(risk_stratum, baseline_attempt_rate)],
    by = "risk_stratum"
  )
  
  # Calculate adjusted rates
  patient_data[, adjusted_rate := ifelse(
    predicted_high_risk,
    baseline_attempt_rate * intervention_rr,
    baseline_attempt_rate
  )]
  
  # Pre-calculate age-specific parameters for all patients
  # Create age grid (current age + 0 to n_cycles)
  max_age <- ceiling(max(patient_data$age)) + n_cycles_lifetime
  min_age <- floor(min(patient_data$age))
  
  # Create age lookup table (one-time calculation)
  age_lookup <- data.table(
    age = min_age:max_age,
    age_mortality = get_background_mortality(min_age:max_age)
  )
  
  if (verbose) cat(sprintf("%.1fs\n", difftime(Sys.time(), step_start, units="secs")))
  
  # =============================================================================
  # STEP 2: INITIALIZE STATE TRACKING (matrix form - much faster than data.table)
  # =============================================================================
  
  if (verbose) cat("  [2/5] Initializing state matrices... ")
  step_start <- Sys.time()
  
  # Use matrices for state probabilities (MUCH faster than data.table for this)
  # Dimensions: [patient, cycle, state]
  n_states <- 3
  state_probs <- array(0, dim = c(n_patients, n_cycles_lifetime + 1, n_states))
  
  # Set initial conditions (all in state 1 at cycle 0)
  state_probs[, 1, 1] <- 1
  
  if (verbose) cat(sprintf("%.1fs\n", difftime(Sys.time(), step_start, units="secs")))
  
  # =============================================================================
  # STEP 3: PRE-CALCULATE TRANSITION PROBABILITIES
  # =============================================================================
  
  if (verbose) cat("  [3/5] Pre-calculating transition probabilities... ")
  step_start <- Sys.time()
  
  # Create matrices for transition probabilities
  # Dimensions: [patient, cycle, transition]
  p_matrices <- list(
    p1_stay = matrix(0, n_patients, n_cycles_lifetime),
    p1_to_prior = matrix(0, n_patients, n_cycles_lifetime),
    p1_to_dead = matrix(0, n_patients, n_cycles_lifetime),
    p2_stay = matrix(0, n_patients, n_cycles_lifetime),
    p2_to_dead = matrix(0, n_patients, n_cycles_lifetime)
  )
  
  # Calculate for all patients and cycles at once
  for (cycle in 0:n_cycles_lifetime-1) {
    
    # Current age for each patient
    current_ages <- patient_data$age + (cycle - 1)
    
    # Get age-dependent mortality (vectorized lookup)
    age_mortality <- age_lookup$age_mortality[match(floor(current_ages), age_lookup$age)]
    # age_mortality <- get_background_mortality(current_ages)  # Just all-cause minus nothing
    
    
    # FROM STATE 1 (no_attempts)
    suicide_attempt_prob <- patient_data$adjusted_rate
    suicide_death_prob <- suicide_attempt_prob * clinical_params$death_per_attempt
    other_death_prob <- age_mortality
    
    # Handle scaling
    total_exit <- suicide_attempt_prob + other_death_prob
    scaling <- ifelse(total_exit > 1, 1.0 / total_exit, 1.0)
    
    suicide_attempt_prob <- suicide_attempt_prob * scaling
    suicide_death_prob <- suicide_death_prob * scaling
    other_death_prob <- other_death_prob * scaling
    
    # Store state 1 transitions
    p_matrices$p1_stay[, cycle] <- 1 - suicide_attempt_prob - other_death_prob
    p_matrices$p1_to_prior[, cycle] <- suicide_attempt_prob - suicide_death_prob
    p_matrices$p1_to_dead[, cycle] <- suicide_death_prob + other_death_prob
    
    # FROM STATE 2 (prior_attempt)
    prior_attempt_prob <- patient_data$adjusted_rate * clinical_params$prior_attempt_multiplier
    prior_death_prob <- prior_attempt_prob * clinical_params$death_per_attempt
    
    total_prior_exit <- prior_death_prob + age_mortality
    prior_scaling <- ifelse(total_prior_exit > 1, 0.99 / total_prior_exit, 1.0)
    
    prior_death_prob <- prior_death_prob * prior_scaling
    age_mortality_prior <- age_mortality * prior_scaling
    
    # Store state 2 transitions
    p_matrices$p2_stay[, cycle] <- 1 - prior_death_prob - age_mortality_prior
    p_matrices$p2_to_dead[, cycle] <- prior_death_prob + age_mortality_prior
  }
  
  if (verbose) cat(sprintf("%.1fs\n", difftime(Sys.time(), step_start, units="secs")))
  
  # =============================================================================
  # STEP 4: MAIN SIMULATION LOOP (VECTORIZED MATRIX OPERATIONS)
  # =============================================================================
  
  if (verbose) cat("  [4/5] Running simulation... ")
  step_start <- Sys.time()
  
  # Track outcomes
  total_attempts <- 0
  total_deaths_suicide <- 0
  total_deaths_other <- 0
  total_person_years <- 0
  
  track_by_cycle <- matrix(0, nrow = n_cycles_lifetime, ncol = 7)
  colnames(track_by_cycle) <- c("cycle", "alive_start", "person_years", "attempts", 
                                "suicide_deaths", "other_deaths", "mean_age")
  
  for (cycle in 1:n_cycles_lifetime) {
    
    # Get current state probabilities (all patients)
    prev_no_attempts <- state_probs[, cycle, 1]
    prev_prior_attempt <- state_probs[, cycle, 2]
    prev_dead <- state_probs[, cycle, 3]
    
    # Apply transitions (VECTORIZED MATRIX OPERATIONS)
    state_probs[, cycle + 1, 1] <- prev_no_attempts * p_matrices$p1_stay[, cycle]
    
    state_probs[, cycle + 1, 2] <- prev_no_attempts * p_matrices$p1_to_prior[, cycle] +
      prev_prior_attempt * p_matrices$p2_stay[, cycle]
    
    state_probs[, cycle + 1, 3] <- prev_no_attempts * p_matrices$p1_to_dead[, cycle] +
      prev_prior_attempt * p_matrices$p2_to_dead[, cycle] +
      prev_dead
    
    # Calculate outcomes (vectorized)
    cycle_attempts <- sum(
      prev_no_attempts * patient_data$adjusted_rate +
        prev_prior_attempt * patient_data$adjusted_rate * clinical_params$prior_attempt_multiplier
    )
    
    cycle_deaths_suicide <- sum(
      prev_no_attempts * patient_data$adjusted_rate * clinical_params$death_per_attempt +
        prev_prior_attempt * patient_data$adjusted_rate * clinical_params$prior_attempt_multiplier * 
        clinical_params$death_per_attempt
    )
    
    # Age-dependent mortality for this cycle
    current_ages <- patient_data$age + (cycle - 1)
    age_mortality_vec <- age_lookup$age_mortality[match(floor(current_ages), age_lookup$age)]
    
    cycle_deaths_other <- sum(
      prev_no_attempts * age_mortality_vec +
        prev_prior_attempt * age_mortality_vec
    )
    
    alive_this_cycle <- sum(prev_no_attempts + prev_prior_attempt)
    
    # Track
    total_attempts <- total_attempts + cycle_attempts
    total_deaths_suicide <- total_deaths_suicide + cycle_deaths_suicide
    total_deaths_other <- total_deaths_other + cycle_deaths_other
    total_person_years <- total_person_years + alive_this_cycle
    
    track_by_cycle[cycle, ] <- c(
      cycle, alive_this_cycle, alive_this_cycle, cycle_attempts,
      cycle_deaths_suicide, cycle_deaths_other, mean(current_ages)
    )
  }
  
  if (verbose) cat(sprintf("%.1fs\n", difftime(Sys.time(), step_start, units="secs")))
  
  # =============================================================================
  # STEP 5: CONVERT TO OUTPUT FORMAT
  # =============================================================================
  
  if (verbose) cat("  [5/5] Converting to output format... ")
  step_start <- Sys.time()
  
  # Convert to hesim stateprobs format (only for non-zero probabilities)
  stateprobs_dt <- data.table()
  strategy_id <- which(strategies$strategy_name == strategy_name)
  
  for (cycle in 0:n_cycles_lifetime) {
    for (state in 1:n_states) {
      # Only keep non-zero probabilities
      nonzero_idx <- which(state_probs[, cycle + 1, state] > 1e-10)
      
      if (length(nonzero_idx) > 0) {
        state_dt <- data.table(
          sample = 1,
          strategy_id = strategy_id,
          patient_id = nonzero_idx,
          grp_id = 1,
          state_id = state,
          t = cycle,
          prob = state_probs[nonzero_idx, cycle + 1, state]
        )
        stateprobs_dt <- rbind(stateprobs_dt, state_dt)
      }
    }
  }
  
  setattr(stateprobs_dt, "class", c("stateprobs", "data.table", "data.frame"))
  
  # Convert track_by_cycle to data.table
  track_dt <- as.data.table(track_by_cycle)
  track_dt[, ':='(
    attempt_rate_per_100k = (attempts / person_years) * 100000,
    suicide_death_rate_per_100k = (suicide_deaths / person_years) * 100000,
    other_death_rate_per_100k = (other_deaths / person_years) * 100000
  )]
  
  if (verbose) cat(sprintf("%.1fs\n", difftime(Sys.time(), step_start, units="secs")))
  
  # Calculate summary
  attempt_rate <- (total_attempts / total_person_years) * 100000
  death_rate <- (total_deaths_suicide / total_person_years) * 100000
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) {
    cat(sprintf("\n  Results: %.1f attempts, %.1f deaths per 100K person-years\n",
                attempt_rate, death_rate))
    cat(sprintf("  ⚡ Total time: %.1f seconds\n", elapsed))
  }
  
  return(list(
    strategy_name = strategy_name,
    stateprobs = stateprobs_dt,
    track_by_cycle = track_dt,
    state_probs_array = state_probs,  # Keep for detailed analysis
    elapsed_time = elapsed,
    summary = list(
      total_attempts = total_attempts,
      total_person_years = total_person_years,
      total_deaths = total_deaths_suicide,
      attempt_rate_per_100k = attempt_rate,
      death_rate_per_100k = death_rate,
      total_deaths_other = total_deaths_other,
      total_deaths_suicide = total_deaths_suicide,
      total_deaths_all = total_deaths_suicide + total_deaths_other,
      suicide_death_rate_per_100k = death_rate,
      other_death_rate_per_100k = (total_deaths_other / total_person_years) * 100000,
      total_death_rate_per_100k = ((total_deaths_suicide + total_deaths_other) / total_person_years) * 100000
    )
  ))
}


# =============================================================================
# 3. OPTIMIZED COST AND QALY CALCULATION
# =============================================================================

calculate_lifetime_costs_qalys_optimized <- function(simulation_result) {
  
  strategy_name <- simulation_result$strategy_name
  state_probs_array <- simulation_result$state_probs_array
  
  cat("Calculating costs and QALYs for", strategy_name, "...\n")
  
  start_time <- Sys.time()
  
  # Get intervention cost
  intervention_cost <- intervention_params$annual_cost[[strategy_name]]
  
  # Pre-calculate background costs by age
  patient_ages <- patients_with_pred$age
  
  # Calculate costs and utilities for all patient-cycles
  total_costs <- 0
  total_qalys <- 0
  
  for (cycle in 0:(n_cycles_lifetime - 1)) {
    
    current_ages <- patient_ages + cycle
    
    # Background costs (vectorized)
    bg_costs <- ifelse(current_ages < 45, cost_params$cost_values$bg_medical_18_44,
                       ifelse(current_ages < 65, cost_params$cost_values$bg_medical_45_64,
                              cost_params$cost_values$bg_medical_65plus))
    
    # State-specific costs
    cost_state_1 <- bg_costs + intervention_cost
    cost_state_2 <- bg_costs + intervention_cost
    cost_state_3 <- 0
    
    # Utilities
    util_state_1 <- clinical_params$base_utility
    util_state_2 <- clinical_params$base_utility * 0.95
    util_state_3 <- 0
    
    # Discount factor
    discount_factor <- 1 / (1 + discount_rate)^cycle
    
    # Calculate discounted costs and QALYs (vectorized)
    cycle_costs <- sum(
      state_probs_array[, cycle + 1, 1] * cost_state_1 +
        state_probs_array[, cycle + 1, 2] * cost_state_2 +
        state_probs_array[, cycle + 1, 3] * cost_state_3
    ) * discount_factor
    
    cycle_qalys <- sum(
      state_probs_array[, cycle + 1, 1] * util_state_1 +
        state_probs_array[, cycle + 1, 2] * util_state_2 +
        state_probs_array[, cycle + 1, 3] * util_state_3
    ) * discount_factor
    
    total_costs <- total_costs + cycle_costs
    total_qalys <- total_qalys + cycle_qalys
  }
  
  cost_per_patient <- total_costs / n_patients
  qalys_per_patient <- total_qalys / n_patients
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  cat("- Cost per patient: $", format(round(cost_per_patient), big.mark = ","), "\n")
  cat("- QALYs per patient:", round(qalys_per_patient, 4), "\n")
  cat(sprintf("- Time: %.1f seconds\n", elapsed))
  
  return(list(
    strategy_name = strategy_name,
    total_costs = total_costs,
    total_qalys = total_qalys,
    cost_per_patient = cost_per_patient,
    qalys_per_patient = qalys_per_patient,
    elapsed_time = elapsed
  ))
}

# =============================================================================
# 4. RUN SIMULATIONS
# =============================================================================

cat("\nRunning OPTIMIZED simulations...\n")
cat(rep("=", 70), "\n")

simulation_start_time <- Sys.time()

# Risk prediction parameters
risk_prediction_params <- list(
  specificity = 0.95,
  sensitivity = 0.25
)

# Apply risk prediction function (same as before)
apply_risk_prediction <- function(patients_dt, sensitivity, specificity) {
  patients_pred <- copy(patients_dt)
  patients_pred <- patients_pred[order(risk_stratum)]
  true_high_risk_threshold <- 950
  patients_pred[, true_high_risk := risk_stratum > true_high_risk_threshold]
  
  n_true_high <- sum(patients_pred$true_high_risk)
  n_true_low <- sum(!patients_pred$true_high_risk)
  
  n_predicted_tp <- round(n_true_high * sensitivity)
  n_predicted_fp <- round(n_true_low * (1 - specificity))
  
  patients_pred[, predicted_high_risk := FALSE]
  
  if (n_predicted_tp > 0) {
    high_risk_ids <- patients_pred[true_high_risk == TRUE]$patient_id
    tp_ids <- sample(high_risk_ids, size = n_predicted_tp, replace = FALSE)
    patients_pred[patient_id %in% tp_ids, predicted_high_risk := TRUE]
  }
  
  if (n_predicted_fp > 0) {
    low_risk_ids <- patients_pred[true_high_risk == FALSE]$patient_id
    fp_ids <- sample(low_risk_ids, size = n_predicted_fp, replace = FALSE)
    patients_pred[patient_id %in% fp_ids, predicted_high_risk := TRUE]
  }
  
  return(patients_pred[, .(patient_id, true_high_risk, predicted_high_risk)])
}

prediction_results <- apply_risk_prediction(patients, 
                                            risk_prediction_params$sensitivity,
                                            risk_prediction_params$specificity)
patients_with_pred <- merge(patients, prediction_results, by = "patient_id")

# Run simulations
all_lifetime_results <- list()
all_economic_results <- list()

for (strategy in strategies$strategy_name) {
  cat("\n", rep("=", 70), "\n", sep="")
  cat("STRATEGY:", strategy, "\n")
  cat(rep("=", 70), "\n", sep="")
  
  sim_result <- run_lifetime_simulation_optimized(strategy, patients_with_pred, verbose = TRUE)
  all_lifetime_results[[strategy]] <- sim_result
  
  econ_result <- calculate_lifetime_costs_qalys_optimized(sim_result)
  all_economic_results[[strategy]] <- econ_result
}

quick_prior_check <- function(simulation_result) {
  
  state_probs_array <- simulation_result$state_probs_array
  strategy_name <- simulation_result$strategy_name
  
  # Check dimensions
  cat(sprintf("\n%s array dimensions: %s\n", strategy_name, 
              paste(dim(state_probs_array), collapse=" × ")))
  
  # Average across all patients and cycles
  # Exclude cycle 0 and final cycle for steady-state estimate
  cycles_to_check <- 10:(n_cycles_lifetime-10)  # Middle cycles for steady state
  
  avg_in_no_attempts <- mean(state_probs_array[, cycles_to_check, 1])
  avg_in_prior_attempt <- mean(state_probs_array[, cycles_to_check, 2])
  avg_dead <- mean(state_probs_array[, cycles_to_check, 3])
  
  # Among living population
  living <- avg_in_no_attempts + avg_in_prior_attempt
  
  if (living > 0) {
    frac_prior <- avg_in_prior_attempt / living
    frac_no <- avg_in_no_attempts / living
    
    cat(sprintf("  Average state distribution (cycles 10-50):\n"))
    cat(sprintf("    No attempts: %.1f%% of living\n", frac_no * 100))
    cat(sprintf("    Prior attempt: %.1f%% of living\n", frac_prior * 100))
    cat(sprintf("    Dead: %.1f%% of total\n", avg_dead * 100))
    
    return(frac_prior)
  } else {
    cat("  ERROR: No living population found\n")
    return(NA)
  }
}


# Run this after simulations
cat("\n=== QUICK STATE DISTRIBUTION CHECK ===\n")
for (strategy in names(all_lifetime_results)) {
  quick_prior_check(all_lifetime_results[[strategy]])
}
simulation_end_time <- Sys.time()
total_time <- as.numeric(difftime(simulation_end_time, simulation_start_time, units = "mins"))

# =============================================================================
# 5. CEA AND RESULTS (same as before)
# =============================================================================

cat("\n", rep("=", 70), "\n", sep="")
cat("COST-EFFECTIVENESS ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep="")

cea_results <- data.table(
  Strategy = names(all_economic_results),
  Attempts_per_100k = sapply(all_lifetime_results, function(x) round(x$summary$attempt_rate_per_100k, 1)),
  Deaths_per_100k = sapply(all_lifetime_results, function(x) round(x$summary$death_rate_per_100k, 1)),
  Cost_per_Patient = sapply(all_economic_results, function(x) round(x$cost_per_patient, 0)),
  QALYs_per_Patient = sapply(all_economic_results, function(x) round(x$qalys_per_patient, 4)),
  Total_Person_Years = sapply(all_lifetime_results, function(x) x$summary$total_person_years),
  Sim_Time_Sec = sapply(all_lifetime_results, function(x) round(x$elapsed_time, 1))
)

print(cea_results)

# Calculate ICERs
baseline_cost <- cea_results[Strategy == "No_Prediction"]$Cost_per_Patient
baseline_qalys <- cea_results[Strategy == "No_Prediction"]$QALYs_per_Patient

cea_results[, Incremental_Cost := Cost_per_Patient - baseline_cost]
cea_results[, Incremental_QALYs := QALYs_per_Patient - baseline_qalys]
cea_results[, ICER := ifelse(Incremental_QALYs > 0, Incremental_Cost / Incremental_QALYs, NA)]

threshold_icer <- 150000
cea_results[, Cost_Effective := ifelse(is.na(ICER), "Baseline", 
                                       ifelse(ICER <= threshold_icer, "Yes", "No"))]

cat("\n=== INCREMENTAL RESULTS ===\n")
print(cea_results[, .(Strategy, Incremental_Cost, Incremental_QALYs, ICER, Cost_Effective)])

cat(sprintf("\n⚡ TOTAL SIMULATION TIME: %.2f minutes\n", total_time))
cat(sprintf("   Average per strategy: %.1f seconds\n", mean(cea_results$Sim_Time_Sec)))

# Save results
if (!dir.exists("output/results")) dir.create("output/results", recursive = TRUE)

save(
  all_lifetime_results, all_economic_results, cea_results,
  total_time, risk_prediction_params, patients_with_pred,
  strategies, patients, states, risk_strata,
  clinical_params, intervention_params,
  n_cycles_lifetime, discount_rate,
  file = "output/results/optimized_lifetime_cea_results_100k_top5.RData"
)

cat("\n✓ Optimized simulation complete!\n")
cat("Results saved to: output/results/optimized_lifetime_cea_results.RData\n")