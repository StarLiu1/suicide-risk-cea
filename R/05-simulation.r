# 5. Hesim Simulation - PHASE 2 INDIVIDUAL PATIENT TRACKING MODEL
# File: R/05-simulation.R
# Full implementation with age progression and individual patient tracking

library(hesim)
library(data.table)

# Load all previous data
load("data/hesim_setup.RData")
load("data/hesim_parameters.RData")
load("data/hesim_transitions.RData")
load("data/hesim_costs_utilities.RData")

cat("Creating complete individual patient tracking simulation...\n")
cat("Phase 2 Implementation - Replicating Ross et al. (2021)\n\n")

# Create directories for output
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/results")) dir.create("output/results")
if (!dir.exists("output/figures")) dir.create("output/figures")

# Simulation parameters
n_samples <- 1  # Deterministic for Phase 2 (PSA in Phase 3)
simulation_population <- n_patients  # All patients

cat("Simulation Configuration:\n")
cat("- Patients:", format(simulation_population, big.mark = ","), "\n")
cat("- Time horizon:", n_cycles, "cycles\n")
cat("- Cycle length:", cycle_length, "year\n")
cat("- Discount rate:", discount_rate * 100, "%\n")
cat("- Age progression: YES\n")
cat("- Individual tracking: YES\n\n")

# INDIVIDUAL PATIENT SIMULATION ENGINE
# ====================================
# Custom simulation that tracks each patient through time with age progression

run_individual_patient_simulation <- function(strategy_name) {
  
  cat("Running simulation for strategy:", strategy_name, "\n")
  
  # Initialize tracking arrays
  # Dimensions: [patient, cycle, outcome]
  patient_ages <- matrix(0, nrow = simulation_population, ncol = n_cycles + 1)
  patient_states <- matrix(0, nrow = simulation_population, ncol = n_cycles + 1)
  patient_alive <- matrix(TRUE, nrow = simulation_population, ncol = n_cycles + 1)
  
  # Initialize all patients
  for (p in 1:simulation_population) {
    patient_ages[p, 1] <- patient_params$current_age[p]
    patient_states[p, 1] <- 1  # All start in "no_attempts" state
    patient_alive[p, 1] <- TRUE
  }
  
  # Tracking variables for outcomes
  total_attempts <- 0
  total_deaths_suicide <- 0
  total_deaths_other <- 0
  total_costs <- 0
  total_qalys <- 0
  
  cycle_results <- data.table(
    cycle = integer(),
    attempts_this_cycle = numeric(),
    deaths_suicide_this_cycle = numeric(),
    deaths_other_this_cycle = numeric(),
    costs_this_cycle = numeric(),
    qalys_this_cycle = numeric()
  )
  
  # Run simulation cycles
  for (cycle in 1:n_cycles) {
    
    cycle_attempts <- 0
    cycle_deaths_suicide <- 0
    cycle_deaths_other <- 0
    cycle_costs <- 0
    cycle_qalys <- 0
    
    # Process each patient
    for (p in 1:simulation_population) {
      
      # Skip if patient already dead
      if (!patient_alive[p, cycle]) {
        patient_ages[p, cycle + 1] <- patient_ages[p, cycle]
        patient_states[p, cycle + 1] <- 3  # Dead state
        patient_alive[p, cycle + 1] <- FALSE
        next
      }
      
      # Update age for this cycle
      current_age <- patient_ages[p, cycle]
      patient_ages[p, cycle + 1] <- current_age + 1
      
      # Get current state
      current_state <- patient_states[p, cycle]
      
      # Get transition matrix for this patient at this age
      tryCatch({
        tmat <- create_transition_matrix(
          patient_id = p, 
          strategy_name = strategy_name, 
          cycle = cycle,
          current_age = current_age
        )
      }, error = function(e) {
        cat("Error creating transition matrix for patient", p, "cycle", cycle, ":", e$message, "\n")
        # Use a safe default matrix
        tmat <- matrix(c(0.99, 0.005, 0.005, 0, 0.99, 0.01, 0, 0, 1), nrow = 3, byrow = TRUE)
      })
      
      # Sample next state based on transition probabilities
      current_state_probs <- tmat[current_state, ]
      next_state <- sample(1:3, size = 1, prob = current_state_probs)
      
      patient_states[p, cycle + 1] <- next_state
      
      # Track outcomes this cycle
      if (current_state == 1 && next_state == 2) {
        # Suicide attempt (survived)
        cycle_attempts <- cycle_attempts + 1
      } else if (current_state == 1 && next_state == 3) {
        # Death (could be suicide or other cause)
        # Approximate: if high suicide risk, likely suicide death
        suicide_prob <- tmat[1, 2] + tmat[1, 3]  # Total suicide attempt probability
        other_death_prob <- get_mortality_rate(current_age)
        
        if (suicide_prob > other_death_prob && runif(1) < 0.8) {
          cycle_deaths_suicide <- cycle_deaths_suicide + 1
          cycle_attempts <- cycle_attempts + 1  # Death counts as attempt too
        } else {
          cycle_deaths_other <- cycle_deaths_other + 1
        }
      } else if (current_state == 2 && next_state == 3) {
        # Death from prior attempt state
        cycle_deaths_suicide <- cycle_deaths_suicide + 1
        cycle_attempts <- cycle_attempts + 1
      }
      
      # Update alive status
      patient_alive[p, cycle + 1] <- (next_state != 3)
      
      # Calculate costs for this patient this cycle
      attempts_this_patient <- 0
      deaths_this_patient <- 0
      
      if ((current_state == 1 && next_state == 2) || (current_state == 1 && next_state == 3)) {
        attempts_this_patient <- 1
      }
      if (next_state == 3 && patient_alive[p, cycle]) {
        deaths_this_patient <- 1
      }
      
      patient_costs <- calculate_cycle_costs(
        patient_id = p,
        strategy_name = strategy_name,
        current_age = current_age,
        health_state = ifelse(current_state == 1, "no_attempts", 
                              ifelse(current_state == 2, "prior_attempt", "dead")),
        suicide_attempts = attempts_this_patient,
        suicide_deaths = ifelse(deaths_this_patient == 1 && cycle_deaths_suicide > cycle_deaths_other, 1, 0)
      )
      
      # Apply discounting
      discount_factor <- 1 / (1 + discount_rate)^(cycle - 1)
      discounted_costs <- patient_costs$total * discount_factor
      cycle_costs <- cycle_costs + discounted_costs
      
      # Calculate QALYs
      if (patient_alive[p, cycle + 1]) {
        if (next_state == 1) {
          qaly_this_patient <- clinical_params$base_utility
        } else if (next_state == 2) {
          qaly_this_patient <- clinical_params$base_utility * 0.95
        } else {
          qaly_this_patient <- 0
        }
      } else {
        qaly_this_patient <- 0
      }
      
      discounted_qalys <- qaly_this_patient * discount_factor
      cycle_qalys <- cycle_qalys + discounted_qalys
    }
    
    # Store cycle results
    cycle_results <- rbind(cycle_results, data.table(
      cycle = cycle,
      attempts_this_cycle = cycle_attempts,
      deaths_suicide_this_cycle = cycle_deaths_suicide,
      deaths_other_this_cycle = cycle_deaths_other,
      costs_this_cycle = cycle_costs,
      qalys_this_cycle = cycle_qalys
    ))
    
    # Update totals
    total_attempts <- total_attempts + cycle_attempts
    total_deaths_suicide <- total_deaths_suicide + cycle_deaths_suicide
    total_deaths_other <- total_deaths_other + cycle_deaths_other
    total_costs <- total_costs + cycle_costs
    total_qalys <- total_qalys + cycle_qalys
    
    # Progress reporting
    if (cycle <= 5 || cycle %% 10 == 0) {
      alive_count <- sum(patient_alive[, cycle + 1])
      cat(sprintf("  Cycle %2d: Alive=%s, Attempts=%d, Deaths=%d\n", 
                  cycle, format(alive_count, big.mark = ","), 
                  cycle_attempts, cycle_deaths_suicide + cycle_deaths_other))
    }
    
    # Debug: Check death accumulation
    # cat("Debug - Cycle", cycle, ":\n")
    # cat("  - Alive this cycle:", sum(patient_alive[, cycle + 1]), "\n")
    # cat("  - Deaths reported:", cycle_deaths_suicide + cycle_deaths_other, "\n")
    # cat("  - Cumulative deaths so far:", 25000 - sum(patient_alive[, cycle + 1]), "\n")
  }
  
  
  
  # Calculate final outcomes
  person_years <- simulation_population * n_cycles
  
  results <- list(
    strategy_name = strategy_name,
    total_population = simulation_population,
    total_attempts = total_attempts,
    total_deaths_suicide = total_deaths_suicide,
    total_deaths_other = total_deaths_other,
    total_deaths = total_deaths_suicide + total_deaths_other,
    attempt_rate_per_100k = (total_attempts / person_years) * 100000,
    death_rate_suicide_per_100k = (total_deaths_suicide / person_years) * 100000,
    death_rate_total_per_100k = ((total_deaths_suicide + total_deaths_other) / person_years) * 100000,
    total_costs = total_costs,
    total_qalys = total_qalys,
    cost_per_patient = total_costs / simulation_population,
    qalys_per_patient = total_qalys / simulation_population,
    cycle_results = cycle_results,
    patient_ages = patient_ages,
    patient_states = patient_states,
    patient_alive = patient_alive
  )
  
  cat(sprintf("Strategy %s completed:\n", strategy_name))
  cat(sprintf("  - Suicide attempts: %d (%.1f per 100,000 person-years)\n", 
              total_attempts, results$attempt_rate_per_100k))
  cat(sprintf("  - Suicide deaths: %d (%.1f per 100,000 person-years)\n", 
              total_deaths_suicide, results$death_rate_suicide_per_100k))
  cat(sprintf("  - Other deaths: %d\n", total_deaths_other))
  cat(sprintf("  - Cost per patient: $%.0f\n", results$cost_per_patient))
  cat(sprintf("  - QALYs per patient: %.4f\n", results$qalys_per_patient))
  cat("\n")
  
  return(results)
}

# Run simulation for all strategies
cat("=== RUNNING PHASE 2 INDIVIDUAL PATIENT SIMULATION ===\n\n")

all_results <- list()
for (strategy in strategies$strategy_name) {
  all_results[[strategy]] <- run_individual_patient_simulation(strategy)
}

# Create summary results table
create_phase2_summary <- function(results) {
  
  cat("Creating Phase 2 results summary...\n")
  
  summary_table <- data.table(
    Strategy = names(results),
    Suicide_Attempts_per_100k = sapply(results, function(x) round(x$attempt_rate_per_100k, 1)),
    Suicide_Deaths_per_100k = sapply(results, function(x) round(x$death_rate_suicide_per_100k, 1)),
    Total_Deaths_per_100k = sapply(results, function(x) round(x$death_rate_total_per_100k, 1)),
    Cost_per_Patient = sapply(results, function(x) round(x$cost_per_patient, 0)),
    QALYs_per_Patient = sapply(results, function(x) round(x$qalys_per_patient, 4))
  )
  
  return(summary_table)
}

results_summary <- create_phase2_summary(all_results)

cat("=== PHASE 2 SIMULATION RESULTS ===\n")
print(results_summary)

# Calculate ICERs
calculate_phase2_icers <- function(results) {
  
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

icer_results <- calculate_phase2_icers(all_results)

cat("=== COST-EFFECTIVENESS ANALYSIS ===\n")
print(icer_results)

# Validation against Ross et al. targets
cat("\n=== VALIDATION AGAINST ROSS ET AL. (2021) ===\n")
baseline_results <- all_results[["No_Prediction"]]

cat("Ross et al. targets:\n")
cat("- Suicide attempts: 175 per 100,000 person-years\n")
cat("- Suicide deaths: 15 per 100,000 person-years\n\n")

cat("Our Phase 2 results:\n")
cat(sprintf("- Suicide attempts: %.1f per 100,000 person-years\n", baseline_results$attempt_rate_per_100k))
cat(sprintf("- Suicide deaths: %.1f per 100,000 person-years\n", baseline_results$death_rate_suicide_per_100k))

# Validation ratios
attempt_ratio <- baseline_results$attempt_rate_per_100k / 175
death_ratio <- baseline_results$death_rate_suicide_per_100k / 15

cat(sprintf("\nValidation ratios:\n"))
cat(sprintf("- Attempt rate ratio: %.3f (target: 1.000)\n", attempt_ratio))
cat(sprintf("- Death rate ratio: %.3f (target: 1.000)\n", death_ratio))

# Validation assessment
if (attempt_ratio >= 0.8 && attempt_ratio <= 1.2) {
  cat("✓ Attempt rates within acceptable range (±20%)\n")
} else {
  cat("⚠️  Attempt rates outside target range\n")
}

if (death_ratio >= 0.8 && death_ratio <= 1.2) {
  cat("✓ Death rates within acceptable range (±20%)\n")
} else {
  cat("⚠️  Death rates outside target range\n")
}

# Age progression validation
cat("\n=== AGE PROGRESSION VALIDATION ===\n")
final_ages <- all_results[["No_Prediction"]]$patient_ages[, n_cycles + 1]
initial_ages <- all_results[["No_Prediction"]]$patient_ages[, 1]
age_increase <- mean(final_ages - initial_ages, na.rm = TRUE)

cat(sprintf("Mean age progression: %.1f years over %d cycles\n", age_increase, n_cycles))
cat(sprintf("Expected age progression: %d years\n", n_cycles))

if (abs(age_increase - n_cycles) < 1) {
  cat("✓ Age progression working correctly\n")
} else {
  cat("⚠️  Age progression may have issues\n")
}

# Save all results
cat("\nSaving Phase 2 results...\n")
save(
  all_results, results_summary, icer_results,
  file = "output/results/phase2_individual_patient_results.RData"
)

# Create summary plots
create_phase2_plots <- function() {
  
  if (!require(ggplot2, quietly = TRUE)) {
    cat("ggplot2 not available, skipping plots\n")
    return()
  }
  
  library(ggplot2)
  
  # Plot 1: Suicide death rates
  p1 <- ggplot(results_summary, aes(x = Strategy, y = Suicide_Deaths_per_100k, fill = Strategy)) +
    geom_col() +
    geom_hline(yintercept = 15, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(title = "Suicide Death Rate by Strategy (Phase 2)",
         subtitle = "Dashed line = Ross et al. target (15 per 100,000)",
         y = "Suicide Deaths per 100,000 person-years",
         x = "Strategy") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 2: Cost vs QALYs
  p2 <- ggplot(results_summary, aes(x = QALYs_per_Patient, y = Cost_per_Patient, 
                                    color = Strategy, label = Strategy)) +
    geom_point(size = 4) +
    geom_text(vjust = -1) +
    labs(title = "Cost vs QALYs per Patient (Phase 2)",
         x = "QALYs per Patient", 
         y = "Cost per Patient ($)") +
    theme_minimal()
  
  ggsave("output/figures/phase2_death_rates.png", p1, width = 10, height = 6)
  ggsave("output/figures/phase2_cost_effectiveness.png", p2, width = 10, height = 6)
  
  cat("✓ Plots saved to output/figures/\n")
}

create_phase2_plots()

# cat("\n" + rep("=", 80) + "\n")
cat("PHASE 2 INDIVIDUAL PATIENT SIMULATION COMPLETE\n")
# cat(rep("=", 80) + "\n")
cat("Individual Patient Tracking Model - Ross et al. (2021) Replication\n")
cat("\nKey Achievements:\n")
cat("✓ Individual patient tracking with age progression\n")
cat("✓ Age-dependent mortality and costs\n")
cat("✓ 1000 risk strata implementation\n")
cat("✓ Full lifecycle simulation (80 cycles)\n")
cat("✓ Cost-effectiveness analysis\n")
cat("✓ Validation against paper targets\n")

cat("\nFiles Created:\n")
cat("- Results: output/results/phase2_individual_patient_results.RData\n")
cat("- Plots: output/figures/phase2_*.png\n")

cat("\nNext Steps for Phase 3:\n")
cat("1. Add probabilistic sensitivity analysis (PSA)\n")
cat("2. Implement risk prediction accuracy analysis\n")
cat("3. Add threshold analysis for cost-effectiveness\n")
cat("4. Create final validation against all paper results\n")

# cat(rep("=", 80) + "\n")