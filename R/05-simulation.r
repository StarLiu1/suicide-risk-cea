# 5. Hesim Simulation and Economic Model
# File: R/05-hesim-simulation.R

library(hesim)
library(data.table)

# Load all previous data
load("data/hesim_setup.RData")
load("data/hesim_parameters.RData")
load("data/hesim_transitions.RData") 
load("data/hesim_costs_utilities.RData")

cat("Creating complete hesim economic model...\n")

# Create directories for output
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/results")) dir.create("output/results")

# Set number of PSA samples (1 for deterministic Phase 1)
n_samples <- 1

cat("Model setup:\n")
cat("- PSA samples:", n_samples, "\n")
cat("- Time horizon:", n_cycles, "cycles\n")
cat("- Cycle length:", cycle_length, "year\n")
cat("- Discount rate:", discount_rate * 100, "%\n")

# Create the transition model using hesim's CohortDtstmTrans
create_hesim_transition_model <- function() {
  
  cat("\nCreating hesim transition model...\n")
  
  # For hesim CohortDtstmTrans, we need to specify transition probabilities
  # We'll create a custom function that returns our pre-calculated matrices
  
  # Create parameter object for transitions
  # We need to format our transition data for hesim
  
  # Create transition probability tables for each strategy/patient combination
  trans_prob_tables <- list()
  
  for (i in 1:nrow(tdata)) {
    row <- tdata[i]
    
    # Get transition matrix
    tmat <- create_transition_matrix(row$stratum, row$strategy_name)
    
    # Convert to probability table format
    prob_table <- data.table(
      strategy_id = row$strategy_id,
      patient_id = row$patient_id,
      from = rep(1:3, each = 3),
      to = rep(1:3, times = 3),
      prob = as.vector(t(tmat))  # Flatten matrix row-wise
    )
    
    trans_prob_tables[[i]] <- prob_table
  }
  
  # Combine all tables
  all_trans_probs <- rbindlist(trans_prob_tables)
  
  # Remove zero probabilities for efficiency
  all_trans_probs <- all_trans_probs[prob > 0]
  
  cat("Transition probabilities table created with", nrow(all_trans_probs), "non-zero transitions\n")
  
  return(all_trans_probs)
}

# Create hesim transition model - simplified approach
create_simple_cohort_model <- function() {
  
  cat("\nCreating simplified cohort model for Phase 1...\n")
  
  # We'll run the simulation manually using our transition matrices
  # and then format results for hesim-style analysis
  
  results <- list()
  
  for (strat_id in 1:nrow(strategies)) {
    
    strategy_name <- strategies$strategy_name[strat_id]
    cat("Simulating strategy:", strategy_name, "\n")
    
    # Initialize population across risk strata
    # Each stratum contributes patients according to population weights
    total_pop <- 10000  # Larger population for more stable results
    
    # Track state populations over time
    # Dimensions: [stratum, cycle, state]
    state_pops <- array(0, dim = c(n_risk_strata, n_cycles + 1, 3))
    
    # Initialize populations
    for (s in 1:n_risk_strata) {
      stratum_pop <- total_pop * risk_strata$population_weight[s]
      state_pops[s, 1, 1] <- stratum_pop  # All start in "no_attempts"
      state_pops[s, 1, 2] <- 0             # None in "prior_attempt"
      state_pops[s, 1, 3] <- 0             # None "dead"
    }
    
    # Run simulation cycles
    total_attempts <- 0
    total_deaths <- 0
    total_costs <- 0
    total_qalys <- 0
    
    for (cycle in 1:n_cycles) {
      
      cycle_attempts <- 0
      cycle_deaths <- 0
      cycle_costs <- 0
      cycle_qalys <- 0
      
      for (stratum in 1:n_risk_strata) {
        
        # Get transition matrix for this stratum and strategy
        tmat <- create_transition_matrix(stratum, strategy_name)
        
        # Current populations
        current_pop <- state_pops[stratum, cycle, ]
        
        # Apply transitions
        new_pop <- as.numeric(current_pop %*% tmat)
        
        # Store new populations
        state_pops[stratum, cycle + 1, ] <- new_pop
        
        # Calculate outcomes this cycle
        new_attempts_stratum <- (new_pop[2] - current_pop[2]) + (new_pop[3] - current_pop[3])
        new_deaths_stratum <- new_pop[3] - current_pop[3]
        
        cycle_attempts <- cycle_attempts + new_attempts_stratum
        cycle_deaths <- cycle_deaths + new_deaths_stratum
        
        # Calculate costs this cycle
        alive_pop <- new_pop[1] + new_pop[2]  # Living population
        
        # Background medical costs
        bg_costs <- alive_pop * cost_params$bg_medical_45_64
        
        # Intervention costs
        interv_costs <- alive_pop * intervention_params$annual_cost[strategy_name]
        
        # Attempt costs
        attempt_costs <- new_attempts_stratum * cost_params$nonfatal_attempt_medical
        death_costs <- new_deaths_stratum * cost_params$fatal_attempt_medical
        
        # Discount costs
        discount_factor <- 1 / (1 + discount_rate)^cycle
        stratum_costs <- (bg_costs + interv_costs + attempt_costs + death_costs) * discount_factor
        cycle_costs <- cycle_costs + stratum_costs
        
        # Calculate QALYs
        qaly_no_attempts <- new_pop[1] * clinical_params$base_utility
        qaly_prior_attempts <- new_pop[2] * clinical_params$base_utility * 0.95  # Slight reduction
        stratum_qalys <- (qaly_no_attempts + qaly_prior_attempts) * discount_factor
        cycle_qalys <- cycle_qalys + stratum_qalys
      }
      
      total_attempts <- total_attempts + cycle_attempts
      total_deaths <- total_deaths + cycle_deaths
      total_costs <- total_costs + cycle_costs
      total_qalys <- total_qalys + cycle_qalys
      
      # Print progress
      if (cycle <= 5 || cycle %% 5 == 0) {
        cat(sprintf("  Cycle %2d: Deaths=%.2f, Attempts=%.2f, Cumulative deaths=%.2f\n", 
                    cycle, cycle_deaths, cycle_attempts, total_deaths))
      }
    }
    
    # Calculate final outcomes
    person_years <- total_pop * n_cycles
    
    strategy_results <- list(
      strategy_id = strat_id,
      strategy_name = strategy_name,
      total_population = total_pop,
      total_attempts = total_attempts,
      total_deaths = total_deaths,
      attempt_rate_per_100k = (total_attempts / person_years) * 100000,
      death_rate_per_100k = (total_deaths / person_years) * 100000,
      total_costs = total_costs,
      total_qalys = total_qalys,
      cost_per_patient = total_costs / total_pop,
      qalys_per_patient = total_qalys / total_pop,
      state_populations = state_pops
    )
    
    results[[strategy_name]] <- strategy_results
    
    cat(sprintf("Strategy %s completed:\n", strategy_name))
    cat(sprintf("  - Death rate: %.2f per 100,000 person-years\n", strategy_results$death_rate_per_100k))
    cat(sprintf("  - Attempt rate: %.2f per 100,000 person-years\n", strategy_results$attempt_rate_per_100k))
    cat(sprintf("  - Cost per patient: $%.0f\n", strategy_results$cost_per_patient))
    cat(sprintf("  - QALYs per patient: %.4f\n\n", strategy_results$qalys_per_patient))
  }
  
  return(results)
}

# Run the simulation
cat("=== RUNNING HESIM PHASE 1 SIMULATION ===\n\n")
simulation_results <- create_simple_cohort_model()

# Create summary table
create_results_summary <- function(results) {
  
  cat("Creating results summary...\n")
  
  summary_table <- data.table(
    Strategy = names(results),
    Death_Rate_per_100k = sapply(results, function(x) round(x$death_rate_per_100k, 2)),
    Attempt_Rate_per_100k = sapply(results, function(x) round(x$attempt_rate_per_100k, 2)),
    Cost_per_Patient = sapply(results, function(x) round(x$cost_per_patient, 0)),
    QALYs_per_Patient = sapply(results, function(x) round(x$qalys_per_patient, 4))
  )
  
  return(summary_table)
}

results_summary <- create_results_summary(simulation_results)

cat("=== SIMULATION RESULTS SUMMARY ===\n")
print(results_summary)

# Calculate ICERs vs baseline
calculate_icers <- function(results) {
  
  cat("\nCalculating ICERs vs No_Prediction baseline...\n")
  
  baseline <- results[["No_Prediction"]]
  baseline_cost <- baseline$cost_per_patient
  baseline_qalys <- baseline$qalys_per_patient
  
  icer_results <- data.table(
    Strategy = names(results),
    Incremental_Cost = numeric(length(results)),
    Incremental_QALYs = numeric(length(results)),
    ICER = character(length(results))
  )
  
  for (i in seq_along(results)) {
    strategy_name <- names(results)[i]
    strategy_results <- results[[strategy_name]]
    
    inc_cost <- strategy_results$cost_per_patient - baseline_cost
    inc_qalys <- strategy_results$qalys_per_patient - baseline_qalys
    
    icer_results[i, Strategy := strategy_name]
    icer_results[i, Incremental_Cost := round(inc_cost, 0)]
    icer_results[i, Incremental_QALYs := round(inc_qalys, 6)]
    
    if (strategy_name == "No_Prediction") {
      icer_results[i, ICER := "Baseline"]
    } else if (inc_qalys <= 0) {
      icer_results[i, ICER := "Dominated"]
    } else {
      icer_value <- inc_cost / inc_qalys
      icer_results[i, ICER := paste0("$", format(round(icer_value), big.mark = ","))]
    }
  }
  
  return(icer_results)
}

icer_results <- calculate_icers(simulation_results)

cat("=== COST-EFFECTIVENESS ANALYSIS ===\n")
print(icer_results)

# Check against paper targets
cat("\n=== VALIDATION AGAINST PAPER TARGETS ===\n")
cat("Target rates from Ross et al. paper:\n")
cat("- Suicide attempts: 175 per 100,000 person-years\n")
cat("- Suicide deaths: 15 per 100,000 person-years\n\n")

baseline_results <- simulation_results[["No_Prediction"]]
cat("Our baseline results:\n")
cat(sprintf("- Suicide attempts: %.2f per 100,000 person-years\n", baseline_results$attempt_rate_per_100k))
cat(sprintf("- Suicide deaths: %.2f per 100,000 person-years\n", baseline_results$death_rate_per_100k))

# Validation metrics
attempt_ratio <- baseline_results$attempt_rate_per_100k / 175
death_ratio <- baseline_results$death_rate_per_100k / 15

cat(sprintf("\nValidation ratios (should be close to 1.0):\n"))
cat(sprintf("- Attempt rate ratio: %.3f\n", attempt_ratio))
cat(sprintf("- Death rate ratio: %.3f\n", death_ratio))

if (death_ratio < 0.1) {
  cat("\n⚠️  WARNING: Death rates are much lower than expected!\n")
  cat("Consider:\n")
  cat("1. Increasing baseline attempt rates in risk_strata\n")
  cat("2. Checking death_per_attempt probability\n")
  cat("3. Using more risk strata to capture intermediate risks\n")
}

# Save all results
save(
  simulation_results, results_summary, icer_results,
  file = "output/results/hesim_phase1_results.RData"
)

cat("\n✓ Hesim Phase 1 simulation completed!\n")
cat("Results saved to output/results/hesim_phase1_results.RData\n")

# Create simple plots
create_basic_plots <- function() {
  
  if (!require(ggplot2, quietly = TRUE)) {
    cat("ggplot2 not available, skipping plots\n")
    return()
  }
  
  library(ggplot2)
  
  if (!dir.exists("output/figures")) dir.create("output/figures")
  
  # Plot 1: Death rates by strategy
  p1 <- ggplot(results_summary, aes(x = Strategy, y = Death_Rate_per_100k, fill = Strategy)) +
    geom_col() +
    labs(title = "Suicide Death Rate by Strategy",
         y = "Deaths per 100,000 person-years",
         x = "Strategy") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 2: Cost vs QALYs
  p2 <- ggplot(results_summary, aes(x = QALYs_per_Patient, y = Cost_per_Patient, 
                                    color = Strategy, label = Strategy)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.5) +
    labs(title = "Cost vs QALYs per Patient",
         x = "QALYs per Patient", 
         y = "Cost per Patient ($)") +
    theme_minimal()
  
  ggsave("output/figures/hesim_death_rates.png", p1, width = 8, height = 6)
  ggsave("output/figures/hesim_cost_effectiveness.png", p2, width = 8, height = 6)
  
  cat("✓ Plots saved to output/figures/\n")
}

create_basic_plots()

cat("\nHesim Phase 1 implementation complete!\n")
cat("Next steps:\n")
cat("1. Review results and validate against paper\n") 
cat("2. Adjust risk parameters if needed\n")
cat("3. Expand to full 1000-stratum model (Phase 2)\n")
cat("4. Add probabilistic sensitivity analysis\n")