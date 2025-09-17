library(glue)

source("R/01-model-setup.R")
source("R/02-parameters.R")
source("R/03-transitions.R")  # Updated version

# Run validation checks first
cat("Running validation checks...\n")

cat("=== TRANSITION MATRIX VALIDATION ===\n\n")

# Test 1: Basic functionality with lowest risk stratum
cat("Test 1: Basic transition matrix creation\n")
test_stratum <- risk_strata[1]  # Lowest risk
tmat <- create_transition_matrix(test_stratum, intervention_effect = 1.0)

cat("✓ Matrix created successfully\n")
cat("Row sums:", rowSums(tmat), "\n")
cat("All probabilities non-negative:", all(tmat >= 0), "\n\n")

# Test 2: Test all risk strata
cat("Test 2: All risk strata validation\n")
test_all_transition_matrices()
cat("\n")

# Test 3: Visualize matrices for different risk levels
cat("Test 3: Matrix visualization\n\n")

cat("LOW RISK STRATUM (1st percentile):\n")
visualize_transition_matrix(stratum_id = 1, intervention_effect = 1.0)
cat("\n")

cat("HIGH RISK STRATUM (10th percentile in our simplified model):\n") 
visualize_transition_matrix(stratum_id = 10, intervention_effect = 1.0)
cat("\n")

cat("HIGH RISK STRATUM with CBT intervention (RR = 0.47):\n")
visualize_transition_matrix(stratum_id = 10, intervention_effect = 0.47)
cat("\n")

# Test 4: Verify intervention effects work as expected
cat("Test 4: Intervention effect verification\n")

high_risk_stratum <- risk_strata[10]  # Highest risk stratum
baseline_rate <- high_risk_stratum$baseline_attempt_rate

cat("Baseline attempt rate:", baseline_rate, "\n")

# Test different interventions
for (intervention_name in names(intervention_effects)) {
  
  effect <- intervention_effects[[intervention_name]]
  tmat <- create_transition_matrix(high_risk_stratum, intervention_effect = effect)
  
  # Calculate effective attempt rate from matrix
  # (probability of transitioning from No_Attempts to Prior_Attempt or Dead)
  effective_rate <- tmat[1, 2] + tmat[1, 3]
  expected_rate <- baseline_rate * effect
  
  cat(sprintf("%-20s: Expected rate = %.6f, Matrix rate = %.6f, Match = %s\n",
              intervention_name, expected_rate, effective_rate, 
              abs(expected_rate - effective_rate) < 1e-10))
}
cat("\n")

# Test 5: Edge case validation
cat("Test 5: Edge cases\n")

# Test with very high attempt rate (should still be valid)
extreme_stratum <- data.frame(
  stratum_id = 999,
  baseline_attempt_rate = 0.5,  # 50% annual attempt rate
  proportion_pop = 0.001
)

tryCatch({
  extreme_tmat <- create_transition_matrix(extreme_stratum, intervention_effect = 1.0)
  cat("✓ Extreme high-risk case handled correctly\n")
  cat("Row sums:", rowSums(extreme_tmat), "\n")
}, error = function(e) {
  cat("✗ Error with extreme case:", e$message, "\n")
})

# Test with zero attempt rate
zero_stratum <- data.frame(
  stratum_id = 0,
  baseline_attempt_rate = 0.0,
  proportion_pop = 0.001
)

tryCatch({
  zero_tmat <- create_transition_matrix(zero_stratum, intervention_effect = 1.0)
  cat("✓ Zero risk case handled correctly\n")
  cat("Matrix:\n")
  print(zero_tmat)
}, error = function(e) {
  cat("✗ Error with zero case:", e$message, "\n")
})

cat("\n=== VALIDATION COMPLETE ===\n")
cat("If all tests passed, your transition matrices are ready for simulation!\n")

# 2. Check parameter consistency
if (length(intervention_effects) != length(intervention_costs)) {
  stop("ERROR: Mismatched intervention parameters!")
}

if (nrow(risk_strata) != n_risk_strata) {
  stop("ERROR: Risk strata count mismatch!")
}

cat("✓ Parameters validated\n")

# 3. Check risk distribution sums to 1
total_prop <- sum(risk_strata$proportion_pop)
if (abs(total_prop - 1.0) > 1e-6) {
  stop(paste("ERROR: Risk strata proportions sum to", total_prop, "not 1.0"))
}

cat("✓ Risk distribution validated\n\n")