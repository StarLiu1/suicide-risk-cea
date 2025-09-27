# 6. Risk Prediction Layer - ADD SENSITIVITY/SPECIFICITY TARGETING
# File: R/06-risk-prediction.R
# Implements risk prediction screening to match Ross et al. Table 2 results

library(data.table)

# Load previous simulation data
load("data/hesim_costs_utilities_100k_1cycle.RData")
# source("R/05-simulation.R")  # Get the run_transition_simulation function
run_transition_simulation <- function(strategy_name, patients_with_pred, verbose = TRUE) {
  
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
        if (runif(1) < intervention_uptake) {  # 99.4% for ACF, 89.9% for CBT
          adjusted_rate <- baseline_rate * intervention_rr
        } else {
          adjusted_rate <- baseline_rate * 1.0  # No intervention effect
        }
      } else {
        # This patient gets no intervention  
        adjusted_rate <- baseline_rate * 1.0
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


cat("=== ADDING RISK PREDICTION LAYER ===\n")
cat("Implementing sensitivity/specificity targeting for interventions\n\n")

# =============================================================================
# 1. RISK PREDICTION PARAMETERS (FROM ROSS ET AL. TABLE 2)
# =============================================================================

# Table 2 uses 95% specificity, 25% sensitivity
risk_prediction_params <- list(
  specificity = 0.95,  # 95% specificity
  sensitivity = 0.25,  # 25% sensitivity
  
  # These will determine who gets interventions
  # Only "predicted high-risk" patients receive interventions
  prediction_threshold_percentile = 0.95  # Top 5% predicted as high-risk (95% specificity)
)

cat("Risk prediction parameters:\n")
cat("- Specificity:", risk_prediction_params$specificity * 100, "%\n")
cat("- Sensitivity:", risk_prediction_params$sensitivity * 100, "%\n\n")

# =============================================================================
# 2. IMPLEMENT RISK PREDICTION ALGORITHM
# =============================================================================

apply_risk_prediction <- function(patients_dt, sensitivity, specificity) {
  
  # Create a copy to avoid modifying original
  patients_pred <- copy(patients_dt)
  
  # Sort patients by true risk (risk stratum)
  patients_pred <- patients_pred[order(risk_stratum)]
  
  # Determine who is truly high-risk vs low-risk
  # Use same threshold as Ross et al: >99th percentile = high-risk
  true_high_risk_threshold <- 990  # >99th percentile (strata 991-1000)
  
  patients_pred[, true_high_risk := risk_stratum > true_high_risk_threshold]
  
  # Count true positives and negatives
  n_true_high <- sum(patients_pred$true_high_risk)
  n_true_low <- sum(!patients_pred$true_high_risk)
  
  cat("True risk distribution:\n")
  cat("- High-risk patients (>99th percentile):", n_true_high, "\n")
  cat("- Low-risk patients (≤99th percentile):", n_true_low, "\n")
  
  # Apply prediction algorithm with specified sensitivity/specificity
  
  # True positives (correctly identified high-risk)
  n_predicted_tp <- round(n_true_high * sensitivity)
  
  # False positives (low-risk incorrectly flagged)  
  n_predicted_fp <- round(n_true_low * (1 - specificity))
  
  # Initialize prediction results
  patients_pred[, predicted_high_risk := FALSE]
  
  # Randomly select true positives from actual high-risk patients
  if (n_predicted_tp > 0) {
    high_risk_ids <- patients_pred[true_high_risk == TRUE]$patient_id
    tp_ids <- sample(high_risk_ids, size = n_predicted_tp, replace = FALSE)
    patients_pred[patient_id %in% tp_ids, predicted_high_risk := TRUE]
  }
  
  # Randomly select false positives from actual low-risk patients
  if (n_predicted_fp > 0) {
    low_risk_ids <- patients_pred[true_high_risk == FALSE]$patient_id
    fp_ids <- sample(low_risk_ids, size = n_predicted_fp, replace = FALSE)
    patients_pred[patient_id %in% fp_ids, predicted_high_risk := TRUE]
  }
  
  # Calculate actual performance metrics
  n_predicted_positive <- sum(patients_pred$predicted_high_risk)
  n_tp_actual <- sum(patients_pred$predicted_high_risk & patients_pred$true_high_risk)
  n_fp_actual <- sum(patients_pred$predicted_high_risk & !patients_pred$true_high_risk)
  
  actual_sensitivity <- n_tp_actual / n_true_high
  actual_specificity <- 1 - (n_fp_actual / n_true_low)
  actual_ppv <- n_tp_actual / n_predicted_positive
  
  cat("\nPrediction performance:\n")
  cat("- Predicted high-risk:", n_predicted_positive, "(", round(n_predicted_positive/nrow(patients_pred)*100, 1), "% of population)\n")
  cat("- Actual sensitivity:", round(actual_sensitivity, 3), "\n")
  cat("- Actual specificity:", round(actual_specificity, 3), "\n")
  cat("- Positive predictive value:", round(actual_ppv, 4), "\n\n")
  
  return(patients_pred[, .(patient_id, risk_stratum, true_high_risk, predicted_high_risk)])
}

# =============================================================================
# 3. MODIFIED SIMULATION WITH RISK PREDICTION
# =============================================================================

run_targeted_simulation <- function(strategy_name, sensitivity, specificity, verbose = TRUE) {
  
  if (verbose) cat("Running targeted simulation:", strategy_name, "\n")
  cat("- Sensitivity:", sensitivity * 100, "%, Specificity:", specificity * 100, "%\n")
  
  # Apply risk prediction to determine who gets intervention
  prediction_results <- apply_risk_prediction(patients, sensitivity, specificity)
  
  # Merge prediction results back to patients
  patients_with_pred <- merge(patients, prediction_results, by = "patient_id")
  
  # For intervention strategies, only apply RR to predicted high-risk patients
  if (strategy_name %in% c("ACF_Intervention", "CBT_Intervention")) {
    
    # Create modified intervention parameters
    base_rr <- intervention_params$rr[[strategy_name]]
    
    # Only predicted high-risk patients get the intervention effect
    n_treated <- sum(patients_with_pred$predicted_high_risk)
    n_total <- nrow(patients_with_pred)
    
    cat("- Patients receiving intervention:", n_treated, "of", n_total, 
        "(", round(n_treated/n_total*100, 1), "%)\n")
    # 
    # # Calculate effective population-level RR
    # # Population RR = (treated_fraction × intervention_RR) + (untreated_fraction × 1.0)
    # treated_fraction <- n_treated / n_total
    # untreated_fraction <- 1 - treated_fraction
    # effective_population_rr <- (treated_fraction * base_rr) + (untreated_fraction * 1.0)
    # 
    # cat("- Base intervention RR:", base_rr, "\n")
    # cat("- Effective population RR:", round(effective_population_rr, 3), "\n")
    # 
    # # Temporarily modify intervention parameters for this simulation
    # original_rr <- intervention_params$rr[[strategy_name]]
    # intervention_params$rr[[strategy_name]] <<- effective_population_rr
    
  } else {
    cat("- No prediction targeting (baseline strategy)\n")
  }
  
  # Run the standard transition simulation
  result <- run_transition_simulation(strategy_name, patients_with_pred, verbose = FALSE)
  
  # Restore original intervention parameters
  # if (strategy_name %in% c("ACF_Intervention", "CBT_Intervention")) {
  #   intervention_params$rr[[strategy_name]] <<- original_rr
  # }
  
  # Add prediction metadata to results
  result$prediction_params <- list(
    sensitivity = sensitivity,
    specificity = specificity,
    n_predicted_positive = sum(prediction_results$predicted_high_risk),
    actual_sensitivity = sum(prediction_results$predicted_high_risk & prediction_results$true_high_risk) / sum(prediction_results$true_high_risk),
    actual_specificity = 1 - sum(prediction_results$predicted_high_risk & !prediction_results$true_high_risk) / sum(!prediction_results$true_high_risk)
  )
    # total_attempts = total_attempts,
    # total_deaths = total_deaths_suicide,
    # attempt_rate_per_100k = attempt_rate,
    # death_rate_per_100k = death_rate,
    # total_deaths_other = total_deaths_other,
    # total_deaths_suicide = total_deaths_suicide,
    # total_deaths_all = total_deaths_suicide + total_deaths_other,
    # # attempt_rate_per_100k = (total_attempts / person_years) * 100000,
    # suicide_death_rate_per_100k = (total_deaths_suicide / person_years) * 100000,
    # other_death_rate_per_100k = (total_deaths_other / person_years) * 100000,
    # total_death_rate_per_100k = ((total_deaths_suicide + total_deaths_other) / person_years) * 100000
  if (verbose) {
    cat("Result:", round(result$summary$attempt_rate_per_100k, 2), "attempts,", 
        round(result$summary$total_attempts, 2), "total attempts \n\n",
        round(result$summary$total_deaths_suicide, 2), "total suicide deaths \n\n",
        round(result$summary$other_death_rate_per_100k, 2), "total other deaths \n\n",
        round(result$summary$suicide_death_rate_per_100k, 2), "suicide deaths per 100K person-years\n\n")
        # round(result$summary$death_rate_per_100k, 2), "suicide deaths per 100K person-years\n\n")
  }
  
  return(result)
}

# =============================================================================
# 4. RUN TARGETED SIMULATIONS (REPLICATE TABLE 2)
# =============================================================================

cat("Running targeted simulations to replicate Ross et al. Table 2...\n\n")

# Parameters from Table 2
table2_sensitivity <- 0.25
table2_specificity <- 0.95

# Run all strategies with risk prediction targeting
targeted_results <- list()

for (strategy in strategies$strategy_name) {
  targeted_results[[strategy]] <- run_targeted_simulation(
    strategy, 
    sensitivity = table2_sensitivity, 
    specificity = table2_specificity,
    verbose = TRUE
  )
}

# =============================================================================
# 5. COMPARE RESULTS TO ROSS ET AL. TABLE 2
# =============================================================================

create_table2_comparison <- function(results) {
  
  cat("=== COMPARISON TO ROSS ET AL. TABLE 2 ===\n\n")
  
  # Ross et al. Table 2 results (baseline and differences)
  ross_baseline <- list(
    attempts = 174.15,
    deaths = 15.34
  )
  
  ross_differences <- list(
    ACF_attempts = -5.90,
    ACF_deaths = -0.52,
    CBT_attempts = -17.76,
    CBT_deaths = -1.56
  )
  
  # Our results
  our_baseline <- results[["No_Prediction"]]$summary
  our_acf <- results[["ACF_Intervention"]]$summary  
  our_cbt <- results[["CBT_Intervention"]]$summary
  
  our_differences <- list(
    ACF_attempts = our_acf$attempt_rate_per_100k - our_baseline$attempt_rate_per_100k,
    ACF_deaths = our_acf$death_rate_per_100k - our_baseline$death_rate_per_100k,
    CBT_attempts = our_cbt$attempt_rate_per_100k - our_baseline$attempt_rate_per_100k,
    CBT_deaths = our_cbt$death_rate_per_100k - our_baseline$death_rate_per_100k
  )
  
  # Create comparison table
  comparison <- data.table(
    Outcome = c("Baseline Attempts", "Baseline Deaths", 
                "ACF Attempt Reduction", "ACF Death Reduction",
                "CBT Attempt Reduction", "CBT Death Reduction"),
    Ross_et_al = c(ross_baseline$attempts, ross_baseline$deaths,
                   ross_differences$ACF_attempts, ross_differences$ACF_deaths,
                   ross_differences$CBT_attempts, ross_differences$CBT_deaths),
    Our_Results = c(our_baseline$attempt_rate_per_100k, our_baseline$death_rate_per_100k,
                    our_differences$ACF_attempts, our_differences$ACF_deaths,
                    our_differences$CBT_attempts, our_differences$CBT_deaths),
    Ratio = c(our_baseline$attempt_rate_per_100k / ross_baseline$attempts,
              our_baseline$death_rate_per_100k / ross_baseline$deaths,
              our_differences$ACF_attempts / ross_differences$ACF_attempts,
              our_differences$ACF_deaths / ross_differences$ACF_deaths,
              our_differences$CBT_attempts / ross_differences$CBT_attempts,
              our_differences$CBT_deaths / ross_differences$CBT_deaths)
  )
  
  comparison[, Ratio := round(Ratio, 3)]
  comparison[, Ross_et_al := round(Ross_et_al, 2)]
  comparison[, Our_Results := round(Our_Results, 2)]
  
  print(comparison)
  
  # Overall assessment
  ratios <- comparison$Ratio
  good_match <- sum(ratios >= 0.8 & ratios <= 1.2, na.rm = TRUE)
  total_outcomes <- length(ratios)
  
  cat("\nValidation Assessment:\n")
  cat("- Outcomes within ±20% of target:", good_match, "of", total_outcomes, "\n")
  
  if (good_match >= 4) {
    cat("✓ GOOD MATCH to Ross et al. Table 2\n")
    validation_status <- "GOOD"
  } else {
    cat("⚠️  PARTIAL MATCH - some outcomes need adjustment\n") 
    validation_status <- "PARTIAL"
  }
  
  return(list(comparison = comparison, validation_status = validation_status))
}

# Create the comparison
table2_comparison <- create_table2_comparison(targeted_results)

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

cat("\nSaving risk prediction results...\n")

save(
  # Risk prediction results
  targeted_results, table2_comparison, risk_prediction_params,
  
  # Functions
  apply_risk_prediction, run_targeted_simulation,
  
  # Keep previous objects
  strategies, patients, states, risk_strata,
  clinical_params, intervention_params,
  
  file = "output/results/risk_prediction_results.RData"
)

cat("\n=== RISK PREDICTION LAYER COMPLETE ===\n")
cat("Key Results:\n")
cat("- Risk prediction targeting implemented\n")
cat("- Table 2 replication:", table2_comparison$validation_status, "\n")
cat("- Sensitivity/specificity effects on intervention targeting working\n")
cat("\nNext: Adjust parameters if needed to improve Table 2 match\n")