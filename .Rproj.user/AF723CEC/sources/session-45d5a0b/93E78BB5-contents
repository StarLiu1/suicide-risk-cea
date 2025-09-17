# 4. Hesim Costs and Utilities
# File: R/04-hesim-costs-utilities.R

library(hesim)
library(data.table)

# Load previous data
load("data/hesim_setup.RData")
load("data/hesim_parameters.RData") 
load("data/hesim_transitions.RData")

cat("Creating hesim cost and utility models...\n")

# Create input data for cost and utility models
# We need data expanded by strategies and states (not transitions)
cost_utility_data <- expand(hesim_dat, by = c("strategies", "patients", "states"))

# Add stratum information
cost_utility_data[, stratum := ceiling(state_id / 3)]  # Each stratum has 3 states

# Add state type information  
cost_utility_data[, state_type := states$state_type[match(state_id, states$state_id)]]

# Add strategy information for intervention costs
cost_utility_data[, intervention_cost := intervention_params$annual_cost[strategy_name]]

cat("Cost/utility input data created with", nrow(cost_utility_data), "rows\n")
cat("Sample of cost/utility data:\n")
print(head(cost_utility_data))

# Create cost model parameters
create_cost_params <- function() {
  
  cat("\nCreating cost model parameters...\n")
  
  # For each row in the data, define costs
  cost_data <- copy(cost_utility_data)
  
  # Background medical costs (state-specific)
  cost_data[state_type == "dead", bg_medical := 0]  # Dead patients have no ongoing costs
  cost_data[state_type != "dead", bg_medical := cost_params$bg_medical_45_64]  # Alive patients
  
  # Intervention costs (only for alive patients)
  cost_data[state_type == "dead", total_intervention_cost := 0]
  cost_data[state_type != "dead", total_intervention_cost := intervention_cost]
  
  # Total annual costs (background + intervention)
  cost_data[, total_annual_cost := bg_medical + total_intervention_cost]
  
  return(cost_data)
}

cost_data <- create_cost_params()

cat("Cost parameters summary:\n")
cat("- Background medical cost (alive):", cost_params$bg_medical_45_64, "\n")
cat("- ACF intervention cost:", intervention_params$annual_cost["ACF_Intervention"], "\n")
cat("- CBT intervention cost:", intervention_params$annual_cost["CBT_Intervention"], "\n")

# Create utility model parameters  
create_utility_params <- function() {
  
  cat("\nCreating utility model parameters...\n")
  
  utility_data <- copy(cost_utility_data)
  
  # Base utilities by state
  utility_data[state_type == "dead", utility := 0]  # Dead = 0 utility
  utility_data[state_type == "no_attempts", utility := clinical_params$base_utility]  # Full utility
  utility_data[state_type == "prior_attempt", utility := clinical_params$base_utility * 0.95]  # Slight reduction
  
  return(utility_data)
}

utility_data <- create_utility_params()

cat("Utility parameters summary:\n")
cat("- No attempts utility:", clinical_params$base_utility, "\n")
cat("- Prior attempt utility:", clinical_params$base_utility * 0.95, "\n")
cat("- Dead utility: 0\n")

# Create hesim cost model using StateVals
create_cost_model <- function(n_samples = 1) {
  
  cat("\nCreating hesim cost model...\n")
  
  # Simple cost model - costs depend on state and strategy
  cost_tbl <- data.table(
    strategy_id = cost_data$strategy_id,
    patient_id = cost_data$patient_id, 
    state_id = cost_data$state_id,
    est = cost_data$total_annual_cost
  )
  
  # Create cost model parameters
  cost_params_obj <- stateval_tbl(
    tbl = cost_tbl,
    dist = "fixed",  # Fixed costs for Phase 1
    hesim_data = hesim_dat
  )
  
  return(cost_params_obj)
}

# Create hesim utility model using StateVals
create_utility_model <- function(n_samples = 1) {
  
  cat("\nCreating hesim utility model...\n")
  
  # Simple utility model - utilities depend on state
  utility_tbl <- data.table(
    strategy_id = utility_data$strategy_id,
    patient_id = utility_data$patient_id,
    state_id = utility_data$state_id, 
    est = utility_data$utility
  )
  
  # Create utility model parameters
  utility_params_obj <- stateval_tbl(
    tbl = utility_tbl,
    dist = "fixed",  # Fixed utilities for Phase 1
    hesim_data = hesim_dat
  )
  
  return(utility_params_obj)
}

# Create the cost and utility parameter objects
cost_params_hesim <- create_cost_model()
utility_params_hesim <- create_utility_model()

cat("✓ Cost and utility models created\n")

# Test the parameter objects
cat("\nTesting parameter objects:\n")
cat("Cost parameters object class:", class(cost_params_hesim), "\n")
cat("Utility parameters object class:", class(utility_params_hesim), "\n")

# Save cost and utility models
save(
  cost_data, utility_data,
  cost_params_hesim, utility_params_hesim,
  file = "data/hesim_costs_utilities.RData"
)

cat("\n✓ Cost and utility models saved to data/hesim_costs_utilities.RData\n")

# Create a function to add transition-specific costs (attempts, deaths)
create_transition_cost_function <- function() {
  
  # This function will be called during simulation to add costs
  # associated with suicide attempts and deaths
  
  transition_cost_function <- function(state_populations, cycle) {
    
    # Calculate new attempts and deaths this cycle
    # This would be called by the main simulation
    
    total_cost <- 0
    
    # Add costs for suicide attempts (non-fatal)
    # Add costs for suicide deaths
    # These would be calculated from changes in state populations
    
    return(total_cost)
  }
  
  return(transition_cost_function)
}

transition_cost_func <- create_transition_cost_function()

save(transition_cost_func, file = "data/hesim_transition_costs.RData")

cat("\nNext: Run 05-hesim-simulation.R to create and run the economic model\n")