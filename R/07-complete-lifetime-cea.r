# Complete Lifetime CEA Simulation - Ross et al. (2021) Replication
# File: R/07-complete-lifetime-cea.R
# Uses existing setup from risk prediction script + proper lifetime tracking

library(hesim)
library(data.table)
library(ggplot2)

# Load the complete setup from risk prediction script
load("data/hesim_costs_utilities.RData")

cat("=== COMPLETE LIFETIME CEA SIMULATION ===\n")
cat("Ross et al. (2021) Suicide Risk Prediction Model\n\n")

# =============================================================================
# 1. SIMULATION CONFIGURATION (LIFETIME HORIZON)
# =============================================================================

# Full lifetime parameters
n_cycles_lifetime <- 1  # ~50 years from mean age 48.8 to end of life
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
  
  intervention_uptake <- intervention_params$uptake[[strategy_name]]
  
  # Run simulation cycles
  for (cycle in 1:n_cycles) {
    
    cycle_attempts <- 0
    cycle_deaths <- 0
    cycle_deaths_other <- 0
    
    # Process each patient
    for (patient in 1:n_patients) {
      
      # Get patient characteristics
      patient_age <- patients$age[patient] + (cycle - 1)  # Age progression
      patient_risk_stratum <- patients$risk_stratum[patient]
      
      # Get baseline attempt rate for this patient
      baseline_rate <- risk_strata[risk_stratum == patient_risk_stratum]$baseline_attempt_rate
      
      # Apply intervention effect
      if (patients_with_pred$predicted_high_risk[patient]) {
        # if (runif(1) < intervention_uptake) {  # 99.4% for ACF, 89.9% for CBT
        #   adjusted_rate <- baseline_rate * intervention_rr 
        # } else {
        #   adjusted_rate <- baseline_rate * 1.0  # No intervention effect
        # }
        adjusted_rate <- baseline_rate * intervention_rr
      } else { 
        # This patient gets no intervention  
        adjusted_rate <- baseline_rate * 1.0
      }
      
      
      # Get age-dependent mortality
      # age_mortality <- get_age_mortality(patient_age)
      age_mortality <- get_background_mortality(patient_age) 
      
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
          scaling <- 1.0 / total_exit
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
          cycle_attempts <- cycle_attempts + current_probs[1] * suicide_death_prob
          cycle_deaths <- cycle_deaths + current_probs[1] * suicide_death_prob
        }
        if (other_death_prob > 0) {
          cycle_deaths_other <- cycle_deaths_other + current_probs[1] * other_death_prob  # OTHER deaths
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
        if (age_mortality_adj > 0) {
          cycle_deaths_other <- cycle_deaths_other + current_probs[2] * age_mortality_adj  # OTHER deaths
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
    total_deaths_other <- total_deaths_other + cycle_deaths_other
    
    # Progress reporting
    if (verbose && (cycle <= 5 || cycle %% 20 == 0)) {
      cat(sprintf("  Cycle %2d: Suicide Deaths=%.2f, Other Deaths=%.2f, Attempts=%.2f\n", 
                  cycle, cycle_deaths, cycle_deaths_other, cycle_attempts))
    }
    
    # Calculate population counts at end of cycle
    alive_count <- sum(stateprobs_array[1, , cycle + 1, 1:2])  # States 1 & 2 (alive)
    dead_count <- sum(stateprobs_array[1, , cycle + 1, 3])     # State 3 (dead)
    prior_attempt_count <- sum(stateprobs_array[1, , cycle + 1, 2])  # State 2 (prior attempts)
    
    # Progress reporting with population counts
    if (verbose && (cycle <= 5 || cycle %% 20 == 0)) {
      cat(sprintf("  Cycle %2d: Alive=%s, Prior Attempts=%s, Dead=%s, New Deaths=%.2f\n", 
                  cycle, 
                  format(round(alive_count), big.mark = ","),
                  format(round(prior_attempt_count), big.mark = ","), 
                  format(round(dead_count), big.mark = ","),
                  cycle_deaths))
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
      death_rate_per_100k = death_rate,
      total_deaths_other = total_deaths_other,
      total_deaths_suicide = total_deaths_suicide,
      total_deaths_all = total_deaths_suicide + total_deaths_other,
      # attempt_rate_per_100k = (total_attempts / person_years) * 100000,
      suicide_death_rate_per_100k = (total_deaths_suicide / person_years) * 100000,
      other_death_rate_per_100k = (total_deaths_other / person_years) * 100000,
      total_death_rate_per_100k = ((total_deaths_suicide + total_deaths_other) / person_years) * 100000
    )
  ))
}

# =============================================================================
# 3. CALCULATE DISCOUNTED COSTS AND QALYS
# =============================================================================

calculate_lifetime_costs_qalys <- function(simulation_result) {
  
  strategy_name <- simulation_result$strategy_name
  stateprobs <- simulation_result$stateprobs
  
  cat("Calculating costs and QALYs for", strategy_name, "...\n")
  
  total_costs <- 0
  total_qalys <- 0
  
  # Cost parameters by strategy
  intervention_cost <- intervention_params$annual_cost[[strategy_name]]
  
  # Process each patient-cycle combination using stateprobs data.table
  for (cycle in 0:(n_cycles_lifetime-1)) {  # Note: cycles 0 to n_cycles-1
    
    discount_factor <- 1 / (1 + discount_rate)^cycle
    
    # Get all state probabilities for this cycle
    cycle_probs <- stateprobs[t == cycle]
    
    if (nrow(cycle_probs) == 0) next  # Skip if no data for this cycle
    
    for (i in 1:nrow(cycle_probs)) {
      
      prob_row <- cycle_probs[i]
      patient_id <- prob_row$patient_id
      state_id <- prob_row$state_id
      prob_value <- prob_row$prob
      
      if (prob_value < 1e-6) next  # Skip negligible probabilities
      
      # Get patient age this cycle
      patient_age <- patients[patient_id == prob_row$patient_id]$age + cycle
      
      # Calculate costs and utilities for this state
      if (state_id == 3) {
        # Dead state - no costs or utilities
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
        if (state_id == 1) {
          state_utility <- clinical_params$base_utility  # No attempts
        } else if (state_id == 2) {
          state_utility <- clinical_params$base_utility * 0.95  # Prior attempt
        } else {
          state_utility <- 0  # Shouldn't happen for living states, but safe default
        }
      }
      
      # Apply discounting and probability weighting
      discounted_cost <- state_cost * discount_factor * prob_value
      discounted_qaly <- state_utility * discount_factor * prob_value
      
      total_costs <- total_costs + discounted_cost
      total_qalys <- total_qalys + discounted_qaly
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