# 5. Hesim Economic Simulation - HYBRID CUSTOM/HESIM APPROACH
# File: R/05-simulation.R
# Custom transition simulation + hesim StateVals + hesim CEA

library(hesim)
library(data.table)

# Load all previous data
load("data/hesim_costs_utilities.RData")

cat("=== HYBRID HESIM/CUSTOM ECONOMIC SIMULATION ===\n")
cat("Combining custom transition logic with hesim StateVals and CEA\n\n")

# Create directories for output
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/results")) dir.create("output/results")
if (!dir.exists("output/figures")) dir.create("output/figures")

# =============================================================================
# 1. SIMULATION CONFIGURATION
# =============================================================================

cat("Simulation configuration:\n")
cat("- Strategies:", nrow(strategies), "\n")
cat("- Patients:", format(n_patients, big.mark = ","), "\n")
cat("- Risk strata:", n_risk_strata, "\n")
cat("- Time horizon:", n_cycles, "cycles\n")
cat("- Cycle length:", cycle_length, "year\n")
cat("- Discount rate:", discount_rate * 100, "%\n\n")

# Simulation parameters
n_samples_psa <- 1  # Deterministic for now (expand for PSA later)
simulation_start_time <- Sys.time()

# =============================================================================
# 2. CUSTOM TRANSITION SIMULATION ENGINE
# =============================================================================

run_transition_simulation <- function(strategy_name, verbose = TRUE) {
  
  if (verbose) cat("Simulating strategy:", strategy_name, "...\n")
  
  # Get intervention parameters
  intervention_rr <- intervention_params$rr[[strategy_name]]
  
  # Initialize state probability tracking
  # Dimensions: [sample, patient, cycle, state]
  n_states <- nrow(states)
  stateprobs_array <- array(0, dim = c(1, n_patients, n_cycles + 1, n_states))
  
  # Initialize all patients in state 1 (no_attempts)
  stateprobs_array[1, , 1, 1] <- 1  # All patients start in state 1
  
  # Track outcomes for validation
  total_attempts <- 0
  total_deaths_suicide <- 0
  total_deaths_other <- 0
  
  # Run simulation cycles
  for (cycle in 1:n_cycles) {
    
    cycle_attempts <- 0
    cycle_deaths <- 0
    
    # Process each patient
    for (patient in 1:n_patients) {
      
      # Get patient characteristics
      patient_age <- patients$age[patient] + (cycle - 1)  # Age progression
      patient_risk_stratum <- patients$risk_stratum[patient]
      
      # Get baseline attempt rate for this patient
      baseline_rate <- risk_strata[risk_stratum == patient_risk_stratum]$baseline_attempt_rate
      
      # Apply intervention effect
      adjusted_rate <- baseline_rate * intervention_rr
      
      # Get age-dependent mortality
      age_mortality <- get_age_mortality(patient_age)
      
      # Current state probabilities for this patient
      current_probs <- stateprobs_array[1, patient, cycle, ]
      
      # Calculate new state probabilities
      new_probs <- rep(0, n_states)
      
      # FROM STATE 1 (no_attempts):
      if (current_probs[1] > 0) {
        
        # Transition probabilities
        suicide_attempt_prob <- adjusted_rate
        suicide_death_prob <- adjusted_rate * clinical_params$death_per_attempt
        other_death_prob <- age_mortality
        
        # Ensure probabilities don't exceed 1
        total_exit <- suicide_attempt_prob + other_death_prob
        if (total_exit > 1) {
          scaling <- 0.99 / total_exit
          suicide_attempt_prob <- suicide_attempt_prob * scaling
          suicide_death_prob <- suicide_death_prob * scaling
          other_death_prob <- other_death_prob * scaling
        }
        
        # Apply transitions
        stay_prob <- 1 - suicide_attempt_prob - other_death_prob
        attempt_survive_prob <- suicide_attempt_prob - suicide_death_prob
        death_prob <- suicide_death_prob + other_death_prob
        
        new_probs[1] <- new_probs[1] + current_probs[1] * stay_prob           # Stay no_attempts
        new_probs[2] <- new_probs[2] + current_probs[1] * attempt_survive_prob # -> prior_attempt
        new_probs[3] <- new_probs[3] + current_probs[1] * death_prob          # -> dead
        
        # Track outcomes
        if (attempt_survive_prob > 0) {
          cycle_attempts <- cycle_attempts + current_probs[1] * attempt_survive_prob
        }
        if (suicide_death_prob > 0) {
          cycle_deaths <- cycle_deaths + current_probs[1] * suicide_death_prob
        }
      }
      
      # FROM STATE 2 (prior_attempt):
      if (current_probs[2] > 0) {
        
        prior_attempt_prob <- adjusted_rate * clinical_params$prior_attempt_multiplier
        prior_death_prob <- prior_attempt_prob * clinical_params$death_per_attempt
        
        total_prior_exit <- prior_death_prob + age_mortality
        if (total_prior_exit > 1) {
          scaling <- 0.99 / total_prior_exit
          prior_death_prob <- prior_death_prob * scaling
          age_mortality_adj <- age_mortality * scaling
        } else {
          age_mortality_adj <- age_mortality
        }
        
        stay_prior_prob <- 1 - prior_death_prob - age_mortality_adj
        death_prior_prob <- prior_death_prob + age_mortality_adj
        
        new_probs[2] <- new_probs[2] + current_probs[2] * stay_prior_prob  # Stay prior_attempt
        new_probs[3] <- new_probs[3] + current_probs[2] * death_prior_prob # -> dead
        
        # Track outcomes
        if (prior_death_prob > 0) {
          cycle_attempts <- cycle_attempts + current_probs[2] * prior_attempt_prob
          cycle_deaths <- cycle_deaths + current_probs[2] * prior_death_prob
        }
      }
      
      # FROM STATE 3 (dead):
      new_probs[3] <- new_probs[3] + current_probs[3]  # Stay dead
      
      # Store new probabilities
      stateprobs_array[1, patient, cycle + 1, ] <- new_probs
    }
    
    # Accumulate outcomes
    total_attempts <- total_attempts + cycle_attempts
    total_deaths_suicide <- total_deaths_suicide + cycle_deaths
    
    # Progress reporting
    if (verbose && (cycle <= 5 || cycle %% 20 == 0)) {
      cat(sprintf("  Cycle %2d: Deaths=%.2f, Attempts=%.2f\n", 
                  cycle, cycle_deaths, cycle_attempts))
    }
  }
  
  # Convert to hesim stateprobs format
  stateprobs_dt <- data.table()
  
  for (sample in 1:1) {
    for (strategy_id in 1:nrow(strategies)) {
      if (strategies$strategy_name[strategy_id] != strategy_name) next
      
      for (patient in 1:n_patients) {
        for (cycle in 0:n_cycles) {
          for (state in 1:n_states) {
            
            prob_val <- stateprobs_array[sample, patient, cycle + 1, state]
            
            if (prob_val > 1e-10) {  # Only store non-zero probabilities
              stateprobs_dt <- rbind(stateprobs_dt, data.table(
                sample = sample,
                strategy_id = strategy_id,
                patient_id = patient,
                grp_id = 1,
                state_id = state,
                t = cycle,
                prob = prob_val
              ))
            }
          }
        }
      }
    }
  }
  
  # Set class for hesim compatibility
  setattr(stateprobs_dt, "class", c("stateprobs", "data.table", "data.frame"))
  
  # Calculate summary statistics
  person_years <- n_patients * n_cycles
  attempt_rate <- (total_attempts / person_years) * 100000
  death_rate <- (total_deaths_suicide / person_years) * 100000
  
  if (verbose) {
    cat(sprintf("  Final results: %.1f attempts, %.1f deaths per 100K person-years\n",
                attempt_rate, death_rate))
  }
  
  return(list(
    strategy_name = strategy_name,
    stateprobs = stateprobs_dt,
    summary = list(
      total_attempts = total_attempts,
      total_deaths = total_deaths_suicide,
      attempt_rate_per_100k = attempt_rate,
      death_rate_per_100k = death_rate
    )
  ))
}

# =============================================================================
# 3. RUN SIMULATIONS FOR ALL STRATEGIES
# =============================================================================

cat("Running transition simulations...\n")

# Store all results
all_simulation_results <- list()
all_stateprobs <- data.table()

# Simulate each strategy
for (i in 1:nrow(strategies)) {
  strategy_name <- strategies$strategy_name[i]
  
  # Run transition simulation
  sim_result <- run_transition_simulation(strategy_name, verbose = TRUE)
  all_simulation_results[[strategy_name]] <- sim_result
  
  # Combine state probabilities
  all_stateprobs <- rbind(all_stateprobs, sim_result$stateprobs)
}

cat("\n✓ All transition simulations completed\n")

# =============================================================================
# 4. CALCULATE COSTS AND QALYS USING HESIM STATEVALS
# =============================================================================

cat("\nCalculating costs and QALYs using hesim StateVals...\n")

# Simulate costs using hesim StateVals
tryCatch({
  costs_result <- cost_model_objects$model$sim(
    stateprobs = all_stateprobs,
    dr = discount_rate
  )
  
  cat("✓ Cost calculations completed\n")
  cat("Cost results dimensions:", dim(costs_result), "\n")
  
}, error = function(e) {
  cat("✗ Error in cost calculations:", e$message, "\n")
  costs_result <- NULL
})

# Simulate utilities using hesim StateVals
tryCatch({
  utilities_result <- utility_model_objects$model$sim(
    stateprobs = all_stateprobs,
    dr = discount_rate
  )
  
  cat("✓ Utility calculations completed\n")
  cat("Utility results dimensions:", dim(utilities_result), "\n")
  
}, error = function(e) {
  cat("✗ Error in utility calculations:", e$message, "\n")
  utilities_result <- NULL
})

# =============================================================================
# 5. CREATE COMPLETE ECONOMIC MODEL (HESIM COHORTDTSTM APPROACH)
# =============================================================================

create_hesim_economic_model <- function() {
  
  cat("\nCreating complete hesim economic model...\n")
  
  # We'll create a hybrid approach:
  # - Use our state probabilities (from custom simulation)
  # - Use hesim StateVals for costs and utilities
  # - Use hesim's summarization and CEA functions
  
  # Create a mock CohortDtstm-style object to hold results
  economic_model <- list(
    stateprobs_ = all_stateprobs,
    costs_ = costs_result,
    qalys_ = utilities_result,
    
    # Add summary method compatible with hesim CEA
    summarize = function(by_grp = FALSE) {
      
      if (is.null(costs_result) || is.null(utilities_result)) {
        cat("Cannot summarize - missing cost or utility results\n")
        return(NULL)
      }
      
      # Calculate mean costs and QALYs by strategy
      cost_summary <- costs_result[, .(costs = sum(costs)), by = .(strategy_id)]
      qaly_summary <- utilities_result[, .(qalys = sum(qalys)), by = .(strategy_id)]
      
      # Merge results
      ce_summary <- merge(cost_summary, qaly_summary, by = "strategy_id")
      
      # Add strategy names
      ce_summary <- merge(ce_summary, strategies[, .(strategy_id, strategy_name)], 
                          by = "strategy_id")
      
      # Set class for hesim CEA compatibility
      setattr(ce_summary, "class", c("ce", "data.table", "data.frame"))
      
      return(ce_summary)
    }
  )
  
  class(economic_model) <- c("custom_hesim_model", "list")
  
  return(economic_model)
}

# Create the economic model
if (!is.null(costs_result) && !is.null(utilities_result)) {
  economic_model <- create_hesim_economic_model()
  
  # Get cost-effectiveness summary
  ce_results <- economic_model$summarize()
  
  cat("\n=== COST-EFFECTIVENESS RESULTS ===\n")
  print(ce_results)
  
} else {
  cat("⚠️  Cannot create economic model due to missing cost/utility results\n")
  economic_model <- NULL
  ce_results <- NULL
}

# =============================================================================
# 6. MANUAL CEA CALCULATIONS (BACKUP APPROACH)
# =============================================================================

create_manual_cea <- function() {
  
  cat("\nCreating manual cost-effectiveness analysis...\n")
  
  # Calculate outcomes by strategy
  cea_results <- data.table(
    Strategy = names(all_simulation_results),
    Attempts_per_100k = sapply(all_simulation_results, function(x) x$summary$attempt_rate_per_100k),
    Deaths_per_100k = sapply(all_simulation_results, function(x) x$summary$death_rate_per_100k),
    Total_Attempts = sapply(all_simulation_results, function(x) x$summary$total_attempts),
    Total_Deaths = sapply(all_simulation_results, function(x) x$summary$total_deaths)
  )
  
  # Add cost and QALY calculations (simplified)
  # We'll use average costs per patient from the StateVals setup
  
  cea_results[, Cost_per_Patient := case_when(
    Strategy == "No_Prediction" ~ mean(cost_model_objects$cost_tbl[strategy_id == 1]$est),
    Strategy == "ACF_Intervention" ~ mean(cost_model_objects$cost_tbl[strategy_id == 2]$est),
    Strategy == "CBT_Intervention" ~ mean(cost_model_objects$cost_tbl[strategy_id == 3]$est),
    TRUE ~ 0
  )]
  
  cea_results[, QALYs_per_Patient := case_when(
    Strategy == "No_Prediction" ~ mean(utility_model_objects$utility_tbl[strategy_id == 1]$est) * n_cycles,
    Strategy == "ACF_Intervention" ~ mean(utility_model_objects$utility_tbl[strategy_id == 2]$est) * n_cycles,
    Strategy == "CBT_Intervention" ~ mean(utility_model_objects$utility_tbl[strategy_id == 3]$est) * n_cycles,
    TRUE ~ 0
  )]
  
  # Calculate ICERs vs baseline
  baseline_cost <- cea_results[Strategy == "No_Prediction"]$Cost_per_Patient
  baseline_qalys <- cea_results[Strategy == "No_Prediction"]$QALYs_per_Patient
  
  cea_results[, Incremental_Cost := Cost_per_Patient - baseline_cost]
  cea_results[, Incremental_QALYs := QALYs_per_Patient - baseline_qalys]
  
  cea_results[, ICER := ifelse(Incremental_QALYs > 0, Incremental_Cost / Incremental_QALYs, NA)]
  
  return(cea_results)
}

# Create manual CEA as backup
manual_cea <- create_manual_cea()

cat("\n=== MANUAL COST-EFFECTIVENESS ANALYSIS ===\n")
print(manual_cea)

# =============================================================================
# 7. VALIDATION AGAINST ROSS ET AL. TARGETS
# =============================================================================

cat("\n=== VALIDATION AGAINST ROSS ET AL. (2021) ===\n")

baseline_results <- all_simulation_results[["No_Prediction"]]$summary

cat("Ross et al. target outcomes:\n")
cat("- Suicide attempts: 175 per 100,000 person-years\n")
cat("- Suicide deaths: 15 per 100,000 person-years\n\n")

cat("Our simulation results:\n")
cat(sprintf("- Suicide attempts: %.1f per 100,000 person-years\n", 
            baseline_results$attempt_rate_per_100k))
cat(sprintf("- Suicide deaths: %.1f per 100,000 person-years\n", 
            baseline_results$death_rate_per_100k))

# Calculate validation ratios
attempt_ratio <- baseline_results$attempt_rate_per_100k / 175
death_ratio <- baseline_results$death_rate_per_100k / 15

cat(sprintf("\nValidation ratios:\n"))
cat(sprintf("- Attempt rate ratio: %.3f (target: 1.000)\n", attempt_ratio))
cat(sprintf("- Death rate ratio: %.3f (target: 1.000)\n", death_ratio))

# Validation assessment
validation_status <- "GOOD"
if (attempt_ratio < 0.5 || attempt_ratio > 2.0) validation_status <- "POOR"
if (death_ratio < 0.5 || death_ratio > 2.0) validation_status <- "POOR"

cat(sprintf("\nValidation status: %s\n", validation_status))

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

cat("\nSaving simulation results...\n")

# Save all results
save(
  # Simulation results
  all_simulation_results, all_stateprobs,
  costs_result, utilities_result,
  
  # Economic analysis
  economic_model, ce_results, manual_cea,
  
  # Validation
  validation_status, attempt_ratio, death_ratio,
  
  # Model objects (for reference)
  cost_model_objects, utility_model_objects, event_cost_functions,
  
  file = "output/results/complete_economic_simulation.RData"
)

# =============================================================================
# 9. CREATE SUMMARY REPORT
# =============================================================================

# cat("\n" + strrep("=", 80) + "\n")
cat("SUICIDE RISK PREDICTION ECONOMIC SIMULATION COMPLETE\n")
cat("Hybrid Custom Transition + Hesim StateVals Approach\n")
# cat(strrep("=", 80) + "\n")

cat("\nModel Configuration:\n")
cat("- Approach: Custom transitions + hesim StateVals + hesim CEA\n")
cat("- Patients:", format(n_patients, big.mark = ","), "\n")
cat("- Risk strata:", n_risk_strata, "\n")
cat("- Time horizon:", n_cycles, "cycles\n")
cat("- Strategies:", paste(strategies$strategy_name, collapse = ", "), "\n")

cat("\nKey Results Summary:\n")
for (strategy in names(all_simulation_results)) {
  result <- all_simulation_results[[strategy]]$summary
  cat(sprintf("- %s: %.1f attempts, %.1f deaths per 100K person-years\n",
              strategy, result$attempt_rate_per_100k, result$death_rate_per_100k))
}

cat("\nValidation vs Ross et al. (2021):", validation_status, "\n")
cat("- Target attempt rate: 175 per 100K (ratio:", round(attempt_ratio, 3), ")\n")
cat("- Target death rate: 15 per 100K (ratio:", round(death_ratio, 3), ")\n")

cat("\nFiles Created:\n")
cat("- Complete results: output/results/complete_economic_simulation.RData\n")

cat("\nNext Steps:\n")
cat("1. Review validation ratios and adjust model if needed\n")
cat("2. Add probabilistic sensitivity analysis (PSA)\n")
cat("3. Create publication-ready figures and tables\n")
cat("4. Compare results to Ross et al. benchmarks\n")

simulation_end_time <- Sys.time()
simulation_duration <- difftime(simulation_end_time, simulation_start_time, units = "mins")
cat(sprintf("\nSimulation completed in %.1f minutes\n", as.numeric(simulation_duration)))

# cat(strrep("=", 80) + "\n")