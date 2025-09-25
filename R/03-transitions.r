# Corrected Hesim Transition Model - Individual Patient Discrete-Time Approach
# This replaces the survival model approach with proper discrete-time transition matrices

library(hesim)
library(data.table)

# Load setup data
load("data/hesim_setup.RData")
load("data/hesim_parameters.RData")

cat("Creating CORRECTED individual patient transition model...\n")
cat("Using discrete-time transition matrices (not survival models)\n\n")

# CORRECTED INDIVIDUAL PATIENT TRANSITION FUNCTION
# =================================================
# This function creates patient-specific transition matrices for each cycle

create_patient_transition_matrix <- function(patient_id, strategy_name, cycle = 1) {
  
  # Ensure we're working with a single patient
  if (length(patient_id) != 1) {
    stop("This function works with one patient at a time. Use vectorized version if needed.")
  }
  patient_params <- merge(patients, 
                      patient_rates[, .(patient_id, baseline_attempt_rate)],
                      by = "patient_id")
  
  # Get patient data (single row)
  patient_row <- patient_params[patient_id == patient_id]
  
  if (nrow(patient_row) == 0) {
    stop(paste("Patient", patient_id, "not found"))
  }
  
  # Extract scalar values (not vectors)
  current_age <- as.numeric(patient_row$age[1]) + (cycle - 1)
  baseline_rate <- as.numeric(patient_row$baseline_attempt_rate[1])
  
  # Get intervention effect (scalar)
  rr <- as.numeric(intervention_params$rr[[strategy_name]])
  
  # Calculate adjusted attempt rate (scalar)
  adjusted_rate <- baseline_rate * rr
  
  # Get age-dependent mortality (scalar)
  age_mortality <- as.numeric(get_age_mortality(current_age))
  
  # Get clinical parameters (scalars)
  death_per_attempt <- as.numeric(clinical_params$death_per_attempt)
  prior_multiplier <- as.numeric(clinical_params$prior_attempt_multiplier)
  
  # Calculate probabilities (all scalars now)
  suicide_attempt_prob <- adjusted_rate
  suicide_death_prob <- adjusted_rate * death_per_attempt
  other_death_prob <- age_mortality
  
  # Ensure probabilities don't exceed 1 (scalar conditional)
  total_exit_prob <- suicide_attempt_prob + other_death_prob
  if (total_exit_prob > 0.99) {
    scale <- 0.99 / total_exit_prob
    suicide_attempt_prob <- suicide_attempt_prob * scale
    suicide_death_prob <- suicide_death_prob * scale  
    other_death_prob <- other_death_prob * scale
  }
  # Create 3x3 transition matrix (now all values are scalars)
  tmat <- matrix(0, nrow = 3, ncol = 3)
  rownames(tmat) <- c("no_attempts", "prior_attempt", "dead")
  colnames(tmat) <- c("no_attempts", "prior_attempt", "dead")
  
  # FROM STATE 1: No attempts
  tmat[1, 1] <- 1 - suicide_attempt_prob - other_death_prob  # Stay no_attempts
  tmat[1, 2] <- suicide_attempt_prob - suicide_death_prob    # Survive attempt -> prior_attempt
  tmat[1, 3] <- suicide_death_prob + other_death_prob        # Die (suicide + other)
  
  # FROM STATE 2: Prior attempt  
  prior_attempt_prob <- adjusted_rate * prior_multiplier
  prior_death_prob <- prior_attempt_prob * death_per_attempt
  
  # Ensure valid probabilities (scalar conditional)
  total_prior_exit <- prior_death_prob + other_death_prob
  if (total_prior_exit > 0.99) {
    scale <- 0.99 / total_prior_exit
    prior_death_prob <- prior_death_prob * scale
    other_death_prob_adj <- other_death_prob * scale
  } else {
    other_death_prob_adj <- other_death_prob
  }
  
  tmat[2, 1] <- 0  # Cannot return to no_attempts
  tmat[2, 2] <- 1 - prior_death_prob - other_death_prob_adj  # Stay prior_attempt
  tmat[2, 3] <- prior_death_prob + other_death_prob_adj      # Die
  
  # FROM STATE 3: Dead (absorbing)
  tmat[3, 1] <- 0
  tmat[3, 2] <- 0
  tmat[3, 3] <- 1
  
  # Validate matrix (scalars only)
  row_sums <- rowSums(tmat)
  if (!all(abs(row_sums - 1) < 1e-8)) {
    cat("Invalid transition matrix for patient", patient_id, "cycle", cycle, "\n")
    cat("Row sums:", row_sums, "\n")
    print(tmat)
    stop("Transition matrix validation failed")
  }
  
  if (any(tmat < 0)) {
    cat("Negative probabilities for patient", patient_id, "cycle", cycle, "\n")
    print(tmat)
    stop("Negative probabilities detected")
  }
  
  return(tmat)
}


# HESIM-COMPATIBLE TRANSITION MODEL CLASS
# ========================================
# Create a custom transition model that works with hesim framework

create_individual_patient_trans_model <- function(hesim_dat, n_samples = 1) {
  
  cat("Creating individual patient transition model for hesim...\n")
  
  # Instead of using hesim's built-in transition models, we'll create a custom approach
  # that can handle individual patient matrices
  
  # Create input data structure
  input_data <- expand(hesim_dat, by = c("strategies", "patients", "states"))
  
  # Add patient parameters to input data
  input_data <- merge(input_data, 
                      patient_rates[, .(patient_id, baseline_attempt_rate)],
                      by = "patient_id")
  
  # Create a custom transition model object
  trans_model <- list(
    input_data = input_data,
    n_samples = n_samples,
    n_states = 3,
    create_matrix_fn = create_patient_transition_matrix,
    
    # Method to simulate state probabilities
    sim_stateprobs = function(n_cycles) {
      
      cat("Simulating state probabilities for", n_cycles, "cycles...\n")
      
      n_patients <- nrow(hesim_dat$patients)
      n_strategies <- nrow(hesim_dat$strategies)
      
      # Initialize state probabilities
      stateprobs <- data.table()
      
      for (sample_id in 1:n_samples) {
        for (strat_id in 1:n_strategies) {
          
          strategy_name <- hesim_dat$strategies$strategy_name[strat_id]
          cat("  Processing strategy:", strategy_name, "\n")
          
          # Initialize patient states (all start in state 1: no_attempts)
          patient_states <- matrix(0, nrow = n_patients, ncol = 3)
          patient_states[, 1] <- 1  # All start in no_attempts state
          
          # Store initial state
          for (state_id in 1:3) {
            initial_prob <- ifelse(state_id == 1, 1, 0)
            
            stateprobs <- rbind(stateprobs, data.table(
              sample = sample_id,
              strategy_id = strat_id,
              patient_id = rep(1:n_patients, each = 1),
              state_id = state_id,
              t = 0,
              prob = rep(initial_prob, n_patients)
            ))
          }
          
          # Simulate through cycles
          for (cycle in 1:n_cycles) {
            
            if (cycle %% 10 == 0) cat("    Cycle", cycle, "\n")
            
            new_patient_states <- matrix(0, nrow = n_patients, ncol = 3)
            
            for (p in 1:n_patients) {
              
              # Get current state probabilities for this patient
              current_state_probs <- patient_states[p, ]
              
              # Get transition matrix for this patient at this cycle
              tmat <- create_patient_transition_matrix(
                patient_id = p,
                strategy_name = strategy_name,
                cycle = cycle
              )
              
              # Apply transition matrix: new_probs = current_probs %*% transition_matrix
              new_state_probs <- current_state_probs %*% tmat
              new_patient_states[p, ] <- as.numeric(new_state_probs)
            }
            
            # Update patient states
            patient_states <- new_patient_states
            
            # Store results for this cycle
            for (state_id in 1:3) {
              stateprobs <- rbind(stateprobs, data.table(
                sample = sample_id,
                strategy_id = strat_id,
                patient_id = 1:n_patients,
                state_id = state_id,
                t = cycle,
                prob = patient_states[, state_id]
              ))
            }
          }
        }
      }
      
      # Store results in the model object
      self$stateprobs_ <- stateprobs
      cat("State probability simulation completed!\n")
      
      return(invisible(self))
    }
  )
  
  # Add self-reference for methods
  environment(trans_model$sim_stateprobs)$self <- trans_model
  
  # Set class for compatibility
  class(trans_model) <- "individual_patient_trans"
  
  return(trans_model)
}

# CREATE THE TRANSITION MODEL
# ============================

cat("Building individual patient transition model...\n")

# Create the transition model
transition_model <- create_individual_patient_trans_model(
  hesim_dat = hesim_dat,
  n_samples = 1  # Deterministic for now
)

cat("Transition model created with", nrow(transition_model$input_data), "parameter combinations\n")

# TEST THE TRANSITION MODEL
# =========================

cat("\nTesting individual patient transition matrices...\n")

# Test a few sample patients
test_patients <- c(1, 100, 1000, 5000)
test_strategies <- c("No_Prediction", "ACF_Intervention", "CBT_Intervention")

for (patient_id in test_patients) {
  if (patient_id <= n_patients) {
    
    patient_params <- merge(patients, 
                            patient_rates[, .(patient_id, baseline_attempt_rate)],
                            by = "patient_id")
    
    
    patient_info <- patient_params[patient_id == patient_id]
    cat(sprintf("\nPatient %d: Age %.1f, Risk stratum %d, Baseline rate %.6f\n",
                patient_id, patient_info$initial_age, patient_info$risk_stratum, 
                patient_info$baseline_attempt_rate))
    
    for (strategy in test_strategies) {
      
      tmat <- create_patient_transition_matrix(
        patient_id = patient_id,
        strategy_name = strategy,
        cycle = 1
      )
      
      attempt_prob <- tmat[1, 2] + tmat[1, 3]  # Total attempt probability
      cat(sprintf("  %s: Attempt prob = %.6f\n", strategy, attempt_prob))
    }
  }
}

# SAVE THE CORRECTED TRANSITION MODEL
# ====================================

cat("\nSaving corrected transition model components...\n")

save(
  create_patient_transition_matrix,
  create_individual_patient_trans_model, 
  transition_model,
  file = "data/corrected_hesim_transitions.RData"
)

# cat("\n" + rep("=", 70) + "\n")
cat("CORRECTED HESIM TRANSITION MODEL COMPLETE\n")
# cat(rep("=", 70) + "\n")
cat("Key Changes Made:\n")
cat("✓ Replaced survival models with discrete-time transition matrices\n")
cat("✓ Individual patient-specific matrices (not averaged)\n") 
cat("✓ Age progression built into matrix calculation\n")
cat("✓ Risk stratum effects preserved per patient\n")
cat("✓ Intervention effects applied individually\n")
cat("✓ Custom hesim-compatible simulation approach\n")

cat("\nNext Steps:\n")
cat("1. Replace your previous transition model with this corrected version\n")
cat("2. Update your simulation script to use transition_model$sim_stateprobs()\n")
cat("3. Test with a small number of cycles first\n")
cat("4. Validate results match expected attempt rates\n")

cat("\nThis approach maintains individual patient heterogeneity while being\n")
cat("compatible with hesim's economic modeling framework.\n")
# cat(rep("=", 70) + "\n")