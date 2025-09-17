# Simulation Diagnostics for Suicide Risk Model

# 1. Check if state probabilities are being tracked correctly
check_state_probabilities <- function(econmod) {
  
  # Get state probabilities from simulation
  state_probs <- econmod$stateprobs_
  
  print("=== STATE PROBABILITY DIAGNOSTICS ===")
  print(paste("Total simulation samples:", length(unique(state_probs$sample))))
  print(paste("Total strategies:", length(unique(state_probs$strategy_id))))
  print(paste("Total states:", length(unique(state_probs$state_id))))
  print(paste("Total cycles:", length(unique(state_probs$t))))
  
  # Check if probabilities sum to 1 for each sample/strategy/time
  prob_sums <- state_probs[, .(prob_sum = sum(prob)), 
                           by = .(sample, strategy_id, t)]
  
  print("Probability sums (should be ~1.0):")
  print(summary(prob_sums$prob_sum))
  
  # Check which states have non-zero probabilities
  nonzero_states <- state_probs[prob > 0.001, 
                                .(avg_prob = mean(prob)), 
                                by = .(state_id)]
  print("States with non-zero probabilities:")
  print(nonzero_states)
  
  return(state_probs)
}

# 2. Check transition patterns specifically for attempts
check_attempt_transitions <- function(econmod) {
  
  print("=== TRANSITION DIAGNOSTICS ===")
  
  # Get the underlying transition model
  trans_model <- econmod$trans_model
  
  # Simulate just the transitions to see what's happening
  if ("CohortDtstmTrans" %in% class(trans_model)) {
    
    # Check a few sample transition matrices
    print("Checking sample transition matrices...")
    
    # This will depend on how hesim stores the transition data
    # Let's check what methods are available
    print("Available methods for transition model:")
    print(methods(class = class(trans_model)))
    
    # Try to access transition probabilities
    if ("sim_stateprobs" %in% methods(class = class(trans_model))) {
      sample_probs <- trans_model$sim_stateprobs(n_cycles = 5)
      print("Sample state probabilities from transition model:")
      print(head(sample_probs, 20))
    }
  }
}

# 3. Check if attempts are being counted correctly
check_attempt_counting <- function(econmod, risk_strata) {
  
  print("=== ATTEMPT COUNTING DIAGNOSTICS ===")
  
  # Get state probabilities
  state_probs <- econmod$stateprobs_
  
  # Identify which states represent "Prior_Attempt" 
  # (these indicate someone made an attempt)
  attempt_states <- states[state_type == "Prior_Attempt"]$state_id
  
  print(paste("Attempt states:", paste(attempt_states, collapse = ", ")))
  
  # Calculate total people in attempt states over time
  attempt_counts <- state_probs[state_id %in% attempt_states, 
                                .(total_in_attempt_states = sum(prob)),
                                by = .(sample, strategy_id, t)]
  
  print("People in attempt states over time (first 10 rows):")
  print(head(attempt_counts, 10))
  
  # Calculate attempt rate per 100k
  total_pop <- state_probs[t == 0, sum(prob)]  # Total population at start
  
  attempt_rate_summary <- attempt_counts[, 
                                         .(avg_attempt_rate_per_100k = mean(total_in_attempt_states) * 100000 / total_pop),
                                         by = .(strategy_id, t)
  ]
  
  print("Attempt rates per 100k by strategy and time:")
  print(head(attempt_rate_summary, 20))
  
  return(list(
    attempt_counts = attempt_counts,
    attempt_rates = attempt_rate_summary
  ))
}

# 4. Manual calculation check
manual_attempt_check <- function(risk_strata) {
  
  print("=== MANUAL CALCULATION CHECK ===")
  
  # Calculate expected attempts based on our risk distribution
  expected_attempts <- risk_strata[, 
                                   .(expected_per_100k = sum(baseline_attempt_rate * proportion_pop * 100000))
  ]
  
  print(paste("Expected attempts per 100k from risk distribution:", 
              expected_attempts$expected_per_100k))
  
  # Check individual strata
  print("Attempt rates by stratum:")
  stratum_rates <- risk_strata[, 
                               .(stratum_id, 
                                 baseline_rate = baseline_attempt_rate,
                                 rate_per_100k = baseline_attempt_rate * 100000,
                                 pop_weight = proportion_pop,
                                 contribution = baseline_attempt_rate * proportion_pop * 100000)
  ]
  print(stratum_rates)
  
  return(stratum_rates)
}

# 5. Check transition matrices for each stratum
check_all_transition_matrices <- function(risk_strata) {
  
  print("=== TRANSITION MATRIX DIAGNOSTICS ===")
  
  for(i in 1:min(3, nrow(risk_strata))) {  # Check first 3 strata
    
    stratum <- risk_strata[i]
    tmat <- create_transition_matrix(stratum)
    
    print(paste("=== Stratum", i, "==="))
    print(paste("Baseline attempt rate:", stratum$baseline_attempt_rate))
    print("Transition matrix:")
    print(round(tmat, 6))
    print(paste("Row sums:", paste(round(rowSums(tmat), 6), collapse = ", ")))
    print(paste("Prob of attempt (1->2):", round(tmat[1,2], 6)))
    print("")
  }
}

# 6. Main diagnostic function
run_full_diagnostics <- function(econmod, risk_strata, states) {
  
  print("##########################################")
  print("RUNNING FULL SIMULATION DIAGNOSTICS")
  print("##########################################")
  
  # Run all diagnostic checks
  state_probs <- check_state_probabilities(econmod)
  check_attempt_transitions(econmod)
  attempt_results <- check_attempt_counting(econmod, risk_strata)
  manual_results <- manual_attempt_check(risk_strata)
  check_all_transition_matrices(risk_strata)
  
  print("##########################################")
  print("DIAGNOSTIC SUMMARY")
  print("##########################################")
  
  return(list(
    state_probs = state_probs,
    attempt_results = attempt_results,
    manual_results = manual_results
  ))
}

# Usage example:
# diagnostic_results <- run_full_diagnostics(econmod, risk_strata, states)