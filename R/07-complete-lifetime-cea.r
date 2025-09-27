# Complete Lifetime CEA Simulation - Ross et al. (2021) Replication
# File: R/07-complete-lifetime-cea.R
# Uses existing setup from risk prediction script + proper lifetime tracking

library(hesim)
library(data.table)
library(ggplot2)

# Load the complete setup from risk prediction script
load("data/hesim_costs_utilities_100k_1cycle.RData")

cat("=== COMPLETE LIFETIME CEA SIMULATION ===\n")
cat("Ross et al. (2021) Suicide Risk Prediction Model\n\n")

# =============================================================================
# 1. SIMULATION CONFIGURATION (LIFETIME HORIZON)
# =============================================================================

# Full lifetime parameters
n_cycles_lifetime <- 50  # ~50 years from mean age 48.8 to end of life
n_samples_psa <- 1       # Deterministic for base case

cat("Lifetime simulation configuration:\n")
cat("- Patients:", format(n_patients, big.mark = ","), "\n")
cat("- Risk strata:", n_risk_strata, "\n")
cat("- Time horizon:", n_cycles_lifetime, "cycles (years)\n")
cat("- Starting age: mean", round(mean(patients$age), 1), "years\n")
cat("- Ending age: mean", round(mean(patients$age) + n_cycles_lifetime, 1), "years\n")
cat("- Discount rate:", discount_rate * 100, "%\n\n")

# =============================================================================
# 2. ENHANCED TRANSITION SIMULATION WITH PROPER PERSON-YEARS TRACKING
# =============================================================================

run_lifetime_simulation <- function(strategy_name, patients_with_pred, verbose = TRUE) {
  
  if (verbose) cat("Running lifetime simulation:", strategy_name, "...\n")
  
  # Get intervention parameters
  intervention_rr <- intervention_params$rr[[strategy_name]]
  intervention_uptake <- intervention_params$uptake[[strategy_name]]
  
  # Initialize tracking arrays
  n_states <- nrow(states)
  stateprobs_array <- array(0, dim = c(1, n_patients, n_cycles_lifetime + 1, n_states))
  
  # Track alive status for proper person-years calculation
  alive_array <- array(TRUE, dim = c(1, n_patients, n_cycles_lifetime + 1))
  
  # Initialize all patients in state 1 (no_attempts) and alive
  stateprobs_array[1, , 1, 1] <- 1  # All start in state 1
  alive_array[1, , 1] <- TRUE       # All start alive
  
  # Enhanced outcome tracking
  lifecycle_outcomes <- data.table(
    cycle = integer(),
    alive_count = numeric(),
    new_attempts = numeric(),
    new_deaths_suicide = numeric(),
    new_deaths_other = numeric(),
    cumulative_attempts = numeric(),
    cumulative_deaths = numeric(),
    person_years_this_cycle = numeric()
  )
  
  total_attempts <- 0
  total_deaths_suicide <- 0
  total_deaths_other <- 0
  total_person_years <- 0
  
  # Run simulation cycles
  for (cycle in 1:n_cycles_lifetime) {
    
    cycle_attempts <- 0
    cycle_deaths_suicide <- 0
    cycle_deaths_other <- 0
    
    # Count alive patients at start of cycle for person-years
    alive_start_cycle <- sum(alive_array[1, , cycle])
    total_person_years <- total_person_years + alive_start_cycle
    
    # Process each patient
    for (patient in 1:n_patients) {
      
      # Skip if patient already dead
      if (!alive_array[1, patient, cycle]) {
        # Dead patients stay in dead state
        stateprobs_array[1, patient, cycle + 1, 3] <- 1
        alive_array[1, patient, cycle + 1] <- FALSE
        next
      }
      
      # Get patient characteristics
      patient_age <- patients$age[patient] + (cycle - 1)  # Age progression
      patient_risk_stratum <- patients$risk_stratum[patient]
      
      # Get baseline attempt rate for this patient
      baseline_rate <- risk_strata[risk_stratum == patient_risk_stratum]$baseline_attempt_rate
      
      # Apply intervention effect if patient is predicted high-risk
      if (patients_with_pred$predicted_high_risk[patient]) {
        if (runif(1) < intervention_uptake) {
          adjusted_rate <- baseline_rate * intervention_rr
        } else {
          adjusted_rate <- baseline_rate  # No intervention due to non-uptake
        }
      } else {
        adjusted_rate <- baseline_rate  # No intervention
      }
      
      # Get age-dependent mortality
      age_mortality <- get_age_mortality(patient_age)
      
      # Current state probabilities for this patient
      current_probs <- stateprobs_array[1, patient, cycle, ]
      
      # Calculate new state probabilities
      new_probs <- rep(0, n_states)
      
      # FROM STATE 1 (no_attempts):
      if (current_probs[1] > 0) {
        
        # Transition probabilities
        suicide_attempt_prob <- pmin(0.8, adjusted_rate)
        suicide_death_prob <- suicide_attempt_prob * clinical_params$death_per_attempt
        other_death_prob <- pmin(0.2, age_mortality)
        
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
        
        # Track outcomes (weight by probability)
        cycle_attempts <- cycle_attempts + current_probs[1] * suicide_attempt_prob
        cycle_deaths_suicide <- cycle_deaths_suicide + current_probs[1] * suicide_death_prob
        cycle_deaths_other <- cycle_deaths_other + current_probs[1] * other_death_prob
      }
      
      # FROM STATE 2 (prior_attempt):
      if (current_probs[2] > 0) {
        
        prior_attempt_prob <- pmin(0.9, adjusted_rate * clinical_params$prior_attempt_multiplier)
        prior_suicide_death_prob <- prior_attempt_prob * clinical_params$death_per_attempt
        prior_other_death_prob <- pmin(0.2, age_mortality)
        
        total_prior_exit <- prior_suicide_death_prob + prior_other_death_prob
        if (total_prior_exit > 1) {
          scaling <- 0.99 / total_prior_exit
          prior_suicide_death_prob <- prior_suicide_death_prob * scaling
          prior_other_death_prob <- prior_other_death_prob * scaling
        }
        
        stay_prior_prob <- 1 - prior_suicide_death_prob - prior_other_death_prob
        death_prior_prob <- prior_suicide_death_prob + prior_other_death_prob
        
        new_probs[2] <- new_probs[2] + current_probs[2] * stay_prior_prob  # Stay prior_attempt
        new_probs[3] <- new_probs[3] + current_probs[2] * death_prior_prob # -> dead
        
        # Track outcomes (weight by probability)
        cycle_attempts <- cycle_attempts + current_probs[2] * prior_attempt_prob
        cycle_deaths_suicide <- cycle_deaths_suicide + current_probs[2] * prior_suicide_death_prob
        cycle_deaths_other <- cycle_deaths_other + current_probs[2] * prior_other_death_prob
      }
      
      # FROM STATE 3 (dead):
      new_probs[3] <- new_probs[3] + current_probs[3]  # Stay dead
      
      # Store new probabilities
      stateprobs_array[1, patient, cycle + 1, ] <- new_probs
      
      # Update alive status
      alive_array[1, patient, cycle + 1] <- (new_probs[3] < 0.5)  # Dead if >50% prob in dead state
    }
    
    # Accumulate total outcomes
    total_attempts <- total_attempts + cycle_attempts
    total_deaths_suicide <- total_deaths_suicide + cycle_deaths_suicide
    total_deaths_other <- total_deaths_other + cycle_deaths_other
    
    # Record cycle outcomes
    alive_end_cycle <- sum(alive_array[1, , cycle + 1])
    
    lifecycle_outcomes <- rbind(lifecycle_outcomes, data.table(
      cycle = cycle,
      alive_count = alive_end_cycle,
      new_attempts = cycle_attempts,
      new_deaths_suicide = cycle_deaths_suicide,
      new_deaths_other = cycle_deaths_other,
      cumulative_attempts = total_attempts,
      cumulative_deaths = total_deaths_suicide + total_deaths_other,
      person_years_this_cycle = alive_start_cycle
    ))
    
    # Progress reporting
    if (verbose && (cycle <= 5 || cycle %% 10 == 0 || cycle == n_cycles_lifetime)) {
      cat(sprintf("  Cycle %2d: Alive=%s, Attempts=%.1f, Deaths=%.1f\n", 
                  cycle, 
                  format(alive_end_cycle, big.mark = ","),
                  cycle_attempts, 
                  cycle_deaths_suicide + cycle_deaths_other))
    }
  }
  
  # Calculate final rates using proper person-years (only from alive patients)
  attempt_rate_per_100k <- (total_attempts / total_person_years) * 100000
  death_rate_per_100k <- (total_deaths_suicide / total_person_years) * 100000
  
  if (verbose) {
    cat(sprintf("  FINAL: %.1f attempts, %.1f deaths per 100K person-years\n",
                attempt_rate_per_100k, death_rate_per_100k))
    cat(sprintf("  Total person-years: %s\n", format(total_person_years, big.mark = ",")))
  }
  
  return(list(
    strategy_name = strategy_name,
    stateprobs = stateprobs_array,
    alive_status = alive_array,
    lifecycle_outcomes = lifecycle_outcomes,
    summary = list(
      total_attempts = total_attempts,
      total_deaths_suicide = total_deaths_suicide,
      total_deaths_other = total_deaths_other,
      total_person_years = total_person_years,
      attempt_rate_per_100k = attempt_rate_per_100k,
      death_rate_per_100k = death_rate_per_100k,
      final_alive = tail(lifecycle_outcomes$alive_count, 1)
    )
  ))
}

# =============================================================================
# 3. CALCULATE DISCOUNTED COSTS AND QALYS
# =============================================================================

calculate_lifetime_costs_qalys <- function(simulation_result) {
  
  strategy_name <- simulation_result$strategy_name
  stateprobs <- simulation_result$stateprobs
  alive_status <- simulation_result$alive_status
  
  cat("Calculating costs and QALYs for", strategy_name, "...\n")
  
  total_costs <- 0
  total_qalys <- 0
  
  # Cost parameters by strategy
  intervention_cost <- intervention_params$annual_cost[[strategy_name]]
  
  # Process each patient-cycle combination
  for (cycle in 1:n_cycles_lifetime) {
    
    discount_factor <- 1 / (1 + discount_rate)^(cycle - 1)
    
    for (patient in 1:n_patients) {
      
      # Skip if patient is dead
      if (!alive_status[1, patient, cycle]) next
      
      # Get patient age this cycle
      patient_age <- patients$age[patient] + (cycle - 1)
      
      # Get state probabilities for this patient this cycle
      state_probs <- stateprobs[1, patient, cycle, ]
      
      # Calculate costs for each state
      for (state in 1:3) {
        
        if (state_probs[state] < 1e-6) next  # Skip negligible probabilities
        
        if (state == 3) {
          # Dead state - no costs
          state_cost <- 0
          state_utility <- 0
        } else {
          # Living states - background + intervention costs
          bg_cost <- case_when(
            patient_age < 45 ~ cost_params$cost_values$bg_medical_18_44,
            patient_age < 65 ~ cost_params$cost_values$bg_medical_45_64,
            TRUE ~ cost_params$cost_values$bg_medical_65plus
          )
          
          state_cost <- bg_cost + intervention_cost
          
          # State-specific utility
          if (state == 1) {
            state_utility <- clinical_params$base_utility  # No attempts
          } else {
            state_utility <- clinical_params$base_utility * 0.95  # Prior attempt
          }
        }
        
        # Add event costs for transitions (simplified)
        # In full model, would track specific transitions and add attempt/death costs
        
        # Apply discounting and probability weighting
        discounted_cost <- state_cost * discount_factor * state_probs[state]
        discounted_qaly <- state_utility * discount_factor * state_probs[state]
        
        total_costs <- total_costs + discounted_cost
        total_qalys <- total_qalys + discounted_qaly
      }
    }
  }
  
  # Calculate per-patient values
  cost_per_patient <- total_costs / n_patients
  qalys_per_patient <- total_qalys / n_patients
  
  cat("- Cost per patient: $", format(round(cost_per_patient), big.mark = ","), "\n")
  cat("- QALYs per patient:", round(qalys_per_patient, 4), "\n")
  
  return(list(
    strategy_name = strategy_name,
    total_costs = total_costs,
    total_qalys = total_qalys,
    cost_per_patient = cost_per_patient,
    qalys_per_patient = qalys_per_patient
  ))
}

# =============================================================================
# 4. RUN COMPLETE LIFETIME SIMULATION FOR ALL STRATEGIES
# =============================================================================

cat("Running complete lifetime simulations...\n\n")

# Risk prediction parameters (from Table 2)
risk_prediction_params <- list(
  specificity = 0.95,
  sensitivity = 0.25
)

# Apply risk prediction (using function from risk prediction script)
apply_risk_prediction <- function(patients_dt, sensitivity, specificity) {
  
  patients_pred <- copy(patients_dt)
  patients_pred <- patients_pred[order(risk_stratum)]
  
  # >99th percentile = high-risk (strata 991-1000)
  true_high_risk_threshold <- 990
  patients_pred[, true_high_risk := risk_stratum > true_high_risk_threshold]
  
  n_true_high <- sum(patients_pred$true_high_risk)
  n_true_low <- sum(!patients_pred$true_high_risk)
  
  # Apply prediction algorithm
  n_predicted_tp <- round(n_true_high * sensitivity)
  n_predicted_fp <- round(n_true_low * (1 - specificity))
  
  patients_pred[, predicted_high_risk := FALSE]
  
  # True positives
  if (n_predicted_tp > 0) {
    high_risk_ids <- patients_pred[true_high_risk == TRUE]$patient_id
    tp_ids <- sample(high_risk_ids, size = n_predicted_tp, replace = FALSE)
    patients_pred[patient_id %in% tp_ids, predicted_high_risk := TRUE]
  }
  
  # False positives
  if (n_predicted_fp > 0) {
    low_risk_ids <- patients_pred[true_high_risk == FALSE]$patient_id
    fp_ids <- sample(low_risk_ids, size = n_predicted_fp, replace = FALSE)
    patients_pred[patient_id %in% fp_ids, predicted_high_risk := TRUE]
  }
  
  return(patients_pred[, .(patient_id, risk_stratum, true_high_risk, predicted_high_risk)])
}

# Apply risk prediction
prediction_results <- apply_risk_prediction(
  patients, 
  risk_prediction_params$sensitivity, 
  risk_prediction_params$specificity
)

# Merge prediction results with patients
patients_with_pred <- merge(patients, prediction_results, by = "patient_id")

# Store all simulation results
all_lifetime_results <- list()
all_economic_results <- list()

# Run simulation for each strategy
for (strategy in strategies$strategy_name) {
  
  cat("=== STRATEGY:", strategy, "===\n")
  
  # Run lifetime simulation
  sim_result <- run_lifetime_simulation(strategy, patients_with_pred, verbose = TRUE)
  all_lifetime_results[[strategy]] <- sim_result
  
  # Calculate costs and QALYs
  econ_result <- calculate_lifetime_costs_qalys(sim_result)
  all_economic_results[[strategy]] <- econ_result
  
  cat("\n")
}

# =============================================================================
# 5. COST-EFFECTIVENESS ANALYSIS
# =============================================================================

cat("=== COST-EFFECTIVENESS ANALYSIS ===\n\n")

# Create CEA results table
cea_results <- data.table(
  Strategy = names(all_economic_results),
  Attempts_per_100k = sapply(all_lifetime_results, function(x) round(x$summary$attempt_rate_per_100k, 1)),
  Deaths_per_100k = sapply(all_lifetime_results, function(x) round(x$summary$death_rate_per_100k, 1)),
  Cost_per_Patient = sapply(all_economic_results, function(x) round(x$cost_per_patient, 0)),
  QALYs_per_Patient = sapply(all_economic_results, function(x) round(x$qalys_per_patient, 4)),
  Total_Person_Years = sapply(all_lifetime_results, function(x) x$summary$total_person_years)
)

print(cea_results)

# Calculate ICERs vs baseline
baseline_cost <- cea_results[Strategy == "No_Prediction"]$Cost_per_Patient
baseline_qalys <- cea_results[Strategy == "No_Prediction"]$QALYs_per_Patient

cea_results[, Incremental_Cost := Cost_per_Patient - baseline_cost]
cea_results[, Incremental_QALYs := QALYs_per_Patient - baseline_qalys]

cea_results[, ICER := ifelse(Incremental_QALYs > 0, 
                             Incremental_Cost / Incremental_QALYs, 
                             NA)]

# Cost-effectiveness assessment at $150,000/QALY threshold
threshold_icer <- 150000
cea_results[, Cost_Effective := ifelse(is.na(ICER), "Baseline", 
                                       ifelse(ICER <= threshold_icer, "Yes", "No"))]

cat("\n=== INCREMENTAL COST-EFFECTIVENESS RESULTS ===\n")
print(cea_results[, .(Strategy, Incremental_Cost, Incremental_QALYs, ICER, Cost_Effective)])

# =============================================================================
# 6. VALIDATION AGAINST ROSS ET AL. TARGETS
# =============================================================================

cat("\n=== VALIDATION AGAINST ROSS ET AL. (2021) ===\n")

baseline_results <- all_lifetime_results[["No_Prediction"]]$summary

cat("Ross et al. targets (Table 1):\n")
cat("- Suicide attempts: 175 per 100,000 person-years\n")
cat("- Suicide deaths: 15 per 100,000 person-years\n\n")

cat("Our lifetime simulation results:\n")
cat(sprintf("- Suicide attempts: %.1f per 100,000 person-years\n", 
            baseline_results$attempt_rate_per_100k))
cat(sprintf("- Suicide deaths: %.1f per 100,000 person-years\n", 
            baseline_results$death_rate_per_100k))

# Validation ratios
attempt_ratio <- baseline_results$attempt_rate_per_100k / 175
death_ratio <- baseline_results$death_rate_per_100k / 15

cat(sprintf("\nValidation ratios (target = 1.000):\n"))
cat(sprintf("- Attempt rate ratio: %.3f\n", attempt_ratio))
cat(sprintf("- Death rate ratio: %.3f\n", death_ratio))

validation_status <- ifelse(
  attempt_ratio >= 0.8 && attempt_ratio <= 1.2 && death_ratio >= 0.8 && death_ratio <= 1.2,
  "EXCELLENT", 
  ifelse(attempt_ratio >= 0.5 && attempt_ratio <= 2.0 && death_ratio >= 0.5 && death_ratio <= 2.0,
         "GOOD", "NEEDS_CALIBRATION")
)

cat("Validation status:", validation_status, "\n")

# =============================================================================
# 7. CREATE SUMMARY VISUALIZATIONS
# =============================================================================

create_cea_plots <- function() {
  
  cat("\nCreating cost-effectiveness visualizations...\n")
  
  # Cost vs QALYs plot
  p1 <- ggplot(cea_results, aes(x = QALYs_per_Patient, y = Cost_per_Patient, 
                                color = Strategy, label = Strategy)) +
    geom_point(size = 4) +
    geom_text(vjust = -1, hjust = 0.5) +
    labs(title = "Cost vs QALYs per Patient (Lifetime Horizon)",
         x = "QALYs per Patient", 
         y = "Cost per Patient ($)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Death rates comparison
  p2 <- ggplot(cea_results, aes(x = Strategy, y = Deaths_per_100k, fill = Strategy)) +
    geom_col() +
    geom_hline(yintercept = 15, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(title = "Suicide Death Rates by Strategy",
         subtitle = "Dashed line = Ross et al. target (15 per 100,000)",
         y = "Suicide Deaths per 100,000 person-years",
         x = "Strategy") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
  
  ggsave("output/figures/lifetime_cost_effectiveness.png", p1, width = 10, height = 8)
  ggsave("output/figures/lifetime_death_rates.png", p2, width = 10, height = 6)
  
  cat("✓ Plots saved to output/figures/\n")
}

create_cea_plots()

# =============================================================================
# 8. SAVE COMPLETE RESULTS
# =============================================================================

cat("\nSaving complete lifetime simulation results...\n")

if (!dir.exists("output/results")) dir.create("output/results", recursive = TRUE)

save(
  # Main results
  all_lifetime_results, all_economic_results, cea_results,
  
  # Validation
  validation_status, attempt_ratio, death_ratio,
  
  # Parameters
  risk_prediction_params, patients_with_pred,
  
  # Model setup (for reference)
  strategies, patients, states, risk_strata,
  clinical_params, intervention_params,
  n_cycles_lifetime, discount_rate,
  
  file = "output/results/complete_lifetime_cea_results.RData"
)

# =============================================================================
# 9. FINAL SUMMARY REPORT
# =============================================================================

cat("\n")
cat("===============================================================================\n")
cat("COMPLETE LIFETIME CEA SIMULATION - ROSS ET AL. (2021) REPLICATION\n")
cat("===============================================================================\n")

cat("\nModel Configuration:\n")
cat("- Population:", format(n_patients, big.mark = ","), "individual patients\n")
cat("- Risk strata: 1000 (individual assignment)\n")
cat("- Time horizon:", n_cycles_lifetime, "years (lifetime)\n")
cat("- Starting age: mean", round(mean(patients$age), 1), "years\n")
cat("- Risk prediction: 95% specificity, 25% sensitivity\n")
cat("- Intervention targeting: Predicted high-risk patients only\n")

cat("\nKey Findings:\n")
for (i in 1:nrow(cea_results)) {
  result <- cea_results[i]
  cat(sprintf("- %s: %.1f attempts, %.1f deaths per 100K person-years\n",
              result$Strategy, result$Attempts_per_100k, result$Deaths_per_100k))
}

cat("\nCost-Effectiveness (vs No Prediction):\n")
for (strategy in cea_results[Strategy != "No_Prediction"]$Strategy) {
  result <- cea_results[Strategy == strategy]
  cost_effective <- result$Cost_Effective
  
  if (!is.na(result$ICER)) {
    cat(sprintf("- %s: ICER = $%s/QALY (%s at $150K threshold)\n",
                strategy, 
                format(round(result$ICER), big.mark = ","),
                cost_effective))
  }
}

cat("\nValidation vs Ross et al. (2021):", validation_status, "\n")
cat("- Attempt rate calibration: ratio =", round(attempt_ratio, 3), "\n")
cat("- Death rate calibration: ratio =", round(death_ratio, 3), "\n")

if (validation_status %in% c("EXCELLENT", "GOOD")) {
  cat("\n✓ Model successfully replicates Ross et al. base case outcomes\n")
  cat("✓ Cost-effectiveness analysis complete\n")
  cat("✓ Ready for sensitivity analysis and policy recommendations\n")
} else {
  cat("\n⚠️  Model calibration needs adjustment to better match targets\n")
  cat("- Consider adjusting risk distribution parameters\n")
  cat("- Review transition probability calculations\n")
}

cat("\nFiles Created:\n")
cat("- Complete results: output/results/complete_lifetime_cea_results.RData\n")
cat("- Figures: output/figures/lifetime_*.png\n")

cat("\nNext Steps:\n")
cat("1. Probabilistic sensitivity analysis (PSA)\n") 
cat("2. Threshold analysis for cost-effectiveness\n")
cat("3. Budget impact analysis\n")
cat("4. Validation against other published benchmarks\n")

cat("\n===============================================================================\n")
cat("SIMULATION COMPLETED SUCCESSFULLY\n")
cat("===============================================================================\n")