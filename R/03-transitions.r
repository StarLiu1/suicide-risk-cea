# 3. Hesim Transition Model (CORRECTED VERSION)
# File: R/03-hesim-transitions.R

library(hesim)
library(data.table)

# Load previous data
load("data/hesim_setup.RData")
load("data/hesim_parameters.RData")

cat("Creating hesim transition model (corrected version)...\n")

# CORRECTED transition matrix function with proper data access
create_transition_matrix <- function(stratum_id, strategy_name, cycle = 1) {
  
  # Validate inputs
  if (stratum_id < 1 || stratum_id > nrow(risk_strata)) {
    stop(paste("Invalid stratum_id:", stratum_id, ". Must be between 1 and", nrow(risk_strata)))
  }
  
  if (!strategy_name %in% names(intervention_params$rr)) {
    stop(paste("Invalid strategy_name:", strategy_name, ". Valid options:", 
               paste(names(intervention_params$rr), collapse = ", ")))
  }
  
  # Get baseline attempt rate using correct data.table indexing
  baseline_rate <- risk_strata$baseline_attempt_rate[stratum_id]
  
  # Get intervention effect using correct list access
  rr <- intervention_params$rr[[strategy_name]]
  
  # Apply intervention effect
  adjusted_rate <- baseline_rate * rr
  
  # Get clinical parameters
  death_prob <- clinical_params$death_per_attempt
  prior_multiplier <- clinical_params$prior_attempt_multiplier
  
  # Validate adjusted rate
  if (adjusted_rate > 1) {
    warning(paste("Attempt rate > 1 for stratum", stratum_id, "strategy", strategy_name, "- clamping to 0.99"))
    adjusted_rate <- 0.99
  }
  
  # Create 3x3 transition matrix
  tmat <- matrix(0, nrow = 3, ncol = 3)
  rownames(tmat) <- c("no_attempts", "prior_attempt", "dead")
  colnames(tmat) <- c("no_attempts", "prior_attempt", "dead")
  
  # STATE 1: No attempts
  tmat[1, 1] <- 1 - adjusted_rate                           # Stay in no_attempts
  tmat[1, 2] <- adjusted_rate * (1 - death_prob)           # Attempt, survive -> prior_attempt  
  tmat[1, 3] <- adjusted_rate * death_prob                  # Attempt, die -> dead
  
  # STATE 2: Prior attempt
  prior_attempt_prob <- adjusted_rate * prior_multiplier
  
  # Ensure prior attempt probability doesn't exceed 1
  if (prior_attempt_prob > 1) {
    prior_attempt_prob <- 1
  }
  
  tmat[2, 1] <- 0                                           # Can't go back to no_attempts
  tmat[2, 2] <- 1 - (prior_attempt_prob * death_prob)      # Stay alive
  tmat[2, 3] <- prior_attempt_prob * death_prob            # Die from attempt
  
  # STATE 3: Dead (absorbing)
  tmat[3, 1] <- 0
  tmat[3, 2] <- 0  
  tmat[3, 3] <- 1
  
  # Validate matrix
  row_sums <- rowSums(tmat)
  
  if (!all(abs(row_sums - 1) < 1e-10)) {
    cat("ERROR: Invalid transition matrix for stratum", stratum_id, "strategy", strategy_name, "\n")
    cat("Matrix:\n")
    print(tmat)
    cat("Row sums:", row_sums, "\n")
    cat("Parameters: baseline_rate=", baseline_rate, ", rr=", rr, ", adjusted_rate=", adjusted_rate, "\n")
    stop("Transition matrix validation failed")
  }
  
  if (any(tmat < 0)) {
    stop("Negative probabilities in transition matrix")
  }
  
  return(tmat)
}

# Test the corrected function
cat("Testing corrected transition matrix creation...\n")

# Test with lowest risk stratum
tryCatch({
  test_matrix <- create_transition_matrix(stratum_id = 1, strategy_name = "No_Prediction")
  cat("✓ Test matrix for lowest risk stratum created successfully:\n")
  print(round(test_matrix, 8))
  cat("Row sums:", rowSums(test_matrix), "\n\n")
}, error = function(e) {
  cat("✗ Error creating test matrix:", e$message, "\n")
  stop("Cannot proceed - fix transition matrix function")
})

# Test with highest risk stratum and CBT intervention
tryCatch({
  test_matrix_high <- create_transition_matrix(stratum_id = n_risk_strata, strategy_name = "CBT_Intervention")
  cat("✓ Test matrix for highest risk stratum with CBT:\n")
  print(round(test_matrix_high, 8))
  cat("Row sums:", rowSums(test_matrix_high), "\n\n")
}, error = function(e) {
  cat("✗ Error creating high-risk matrix:", e$message, "\n")
  stop("Cannot proceed - fix transition matrix function")
})

# Create transition list for all combinations (corrected version)
create_transition_list <- function() {
  
  cat("Creating transition matrices for all combinations...\n")
  
  transition_list <- list()
  counter <- 1
  
  for (strat_id in 1:nrow(strategies)) {
    strategy_name <- strategies$strategy_name[strat_id]
    
    for (patient_id in 1:nrow(patients)) {
      # In Phase 1, patient_id maps directly to risk stratum
      stratum_id <- ((patient_id - 1) %% n_risk_strata) + 1
      
      # Create transition matrix for this combination
      tmat <- create_transition_matrix(stratum_id, strategy_name)
      
      # Store with metadata
      transition_list[[counter]] <- list(
        strategy_id = strat_id,
        patient_id = patient_id,
        stratum_id = stratum_id,
        strategy_name = strategy_name,
        transition_matrix = tmat
      )
      
      counter <- counter + 1
    }
  }
  
  cat("✓ Created", length(transition_list), "transition matrices\n")
  return(transition_list)
}

# Create all transition matrices
transition_matrices <- create_transition_list()

# Create input data for transition model (corrected version)
tdata <- expand(hesim_dat, by = c("strategies", "patients"))

# Add stratum mapping
tdata[, stratum := ((patient_id - 1) %% n_risk_strata) + 1]

# Add baseline rates and intervention effects using corrected syntax
tdata[, baseline_attempt_rate := risk_strata$baseline_attempt_rate[stratum]]
tdata[, rr := sapply(strategy_name, function(x) intervention_params$rr[[x]])]
tdata[, adjusted_rate := baseline_attempt_rate * rr]

cat("✓ Transition input data created with", nrow(tdata), "rows\n")

# Display sample of transition data
cat("Sample of transition data:\n")
print(head(tdata, 10))

# Save transition model components
save(
  transition_matrices, tdata, create_transition_matrix,
  file = "data/hesim_transitions.RData"
)

cat("\n✓ Corrected transition model saved to data/hesim_transitions.RData\n")

# Final validation - test matrices give reasonable death rates
cat("\nFinal validation - expected death rates:\n")

for (stratum in c(1, 5, 10)) {
  
  stratum_rate <- risk_strata$baseline_attempt_rate[stratum]
  stratum_pop_weight <- risk_strata$population_weight[stratum]
  
  cat(sprintf("\nStratum %d (rate=%.6f, pop_weight=%.2f):\n", 
              stratum, stratum_rate, stratum_pop_weight))
  
  for (strategy in c("No_Prediction", "ACF_Intervention", "CBT_Intervention")) {
    
    tmat <- create_transition_matrix(stratum, strategy)
    
    # Annual death probability from "no attempts" state
    annual_death_prob <- tmat[1, 3]
    
    # Expected deaths per 100,000 for this stratum
    deaths_per_100k <- annual_death_prob * 100000
    
    cat(sprintf("  %s: %.8f annual death prob (%.2f per 100,000)\n", 
                strategy, annual_death_prob, deaths_per_100k))
  }
}

# Calculate overall expected death rate
total_expected_death_rate <- 0
for (stratum in 1:n_risk_strata) {
  tmat <- create_transition_matrix(stratum, "No_Prediction")
  death_prob <- tmat[1, 3]
  pop_weight <- risk_strata$population_weight[stratum]
  contribution <- death_prob * pop_weight * 100000
  total_expected_death_rate <- total_expected_death_rate + contribution
}

cat(sprintf("\nExpected overall death rate: %.2f per 100,000 person-years\n", total_expected_death_rate))
cat("Target from paper: 15.0 per 100,000 person-years\n")

if (total_expected_death_rate < 5) {
  cat("⚠️  Death rate is low - consider increasing baseline attempt rates\n")
} else if (total_expected_death_rate > 50) {
  cat("⚠️  Death rate is high - consider decreasing baseline attempt rates\n") 
} else {
  cat("✓ Death rate is in reasonable range\n")
}

cat("\nNext: Run 04-hesim-costs-utilities.R to define cost and utility models\n")