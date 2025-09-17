# Master Hesim Runner Script
# File: run_hesim_phase1.R

# This script runs all hesim components in sequence for Phase 1

cat("=== HESIM PHASE 1 SUICIDE RISK PREDICTION MODEL ===\n")
cat("Based on Ross et al. (2021) JAMA Psychiatry\n\n")

# Clear workspace and set up environment
rm(list = ls())

# Load required packages
required_packages <- c("hesim", "data.table", "ggplot2")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

cat("Required packages loaded: hesim, data.table, ggplot2\n\n")

# Create directory structure
directories <- c("R", "data", "output", "output/results", "output/figures")
for (dir in directories) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
}

# Define the sequence of scripts to run
scripts <- c(
  "01-model-setup.r",
  "02-parameters.r", 
  "03-transitions.r",
  "04-costs-utilities.r",
  "05-simulation.r"
)

cat("\n=== RUNNING HESIM COMPONENTS ===\n")

# Run each script in sequence
for (i in seq_along(scripts)) {
  
  script_name <- scripts[i]
  script_path <- file.path("R", script_name)
  
  cat(sprintf("\n--- Step %d: %s ---\n", i, script_name))
  
  # Check if script exists
  if (!file.exists(script_path)) {
    cat("⚠️  Script not found:", script_path, "\n")
    cat("Please make sure you've saved all the hesim component scripts to the R/ directory\n")
    cat("Expected scripts:\n")
    for (s in scripts) {
      cat("  - R/", s, "\n")
    }
    stop("Missing required scripts")
  }
  
  # Run the script
  tryCatch({
    source(script_path)
    cat("✓ Completed:", script_name, "\n")
  }, error = function(e) {
    cat("✗ Error in", script_name, ":", e$message, "\n")
    stop(paste("Failed at step", i, "-", script_name))
  })
}

cat("\n=== HESIM PHASE 1 COMPLETE ===\n")

# Load final results and display summary
if (file.exists("output/results/hesim_phase1_results.RData")) {
  
  load("output/results/hesim_phase1_results.RData")
  
  cat("\n=== FINAL RESULTS SUMMARY ===\n")
  print(results_summary)
  
  cat("\n=== COST-EFFECTIVENESS ANALYSIS ===\n")
  print(icer_results)
  
  # Check if results are reasonable
  baseline <- simulation_results[["No_Prediction"]]
  
  cat("\n=== VALIDATION CHECK ===\n")
  cat("Baseline results vs Ross et al. targets:\n")
  cat(sprintf("Death rate: %.2f vs 15.0 per 100,000 (target)\n", baseline$death_rate_per_100k))
  cat(sprintf("Attempt rate: %.2f vs 175.0 per 100,000 (target)\n", baseline$attempt_rate_per_100k))
  
  death_ratio <- baseline$death_rate_per_100k / 15.0
  attempt_ratio <- baseline$attempt_rate_per_100k / 175.0
  
  if (death_ratio > 0.5 && death_ratio < 2.0) {
    cat("✓ Death rates are in reasonable range\n")
  } else {
    cat("⚠️  Death rates may need calibration\n")
  }
  
  if (attempt_ratio > 0.5 && attempt_ratio < 2.0) {
    cat("✓ Attempt rates are in reasonable range\n")
  } else {
    cat("⚠️  Attempt rates may need calibration\n")
  }
  
  # Check for cost-effectiveness
  acf_icer <- icer_results[Strategy == "ACF_Intervention"]$ICER
  cbt_icer <- icer_results[Strategy == "CBT_Intervention"]$ICER
  
  cat("\nCost-effectiveness at $150,000/QALY threshold:\n")
  
  # Extract numeric ICER values for comparison
  if (acf_icer != "Dominated" && acf_icer != "Baseline") {
    acf_numeric <- as.numeric(gsub("[^0-9]", "", acf_icer))
    if (acf_numeric <= 150000) {
      cat("✓ ACF intervention is cost-effective\n")
    } else {
      cat("✗ ACF intervention exceeds cost-effectiveness threshold\n")
    }
  }
  
  if (cbt_icer != "Dominated" && cbt_icer != "Baseline") {
    cbt_numeric <- as.numeric(gsub("[^0-9]", "", cbt_icer))
    if (cbt_numeric <= 150000) {
      cat("✓ CBT intervention is cost-effective\n")
    } else {
      cat("✗ CBT intervention exceeds cost-effectiveness threshold\n")
    }
  }
  
} else {
  cat("⚠️  Results file not found - simulation may have failed\n")
}

cat("\n=== FILES CREATED ===\n")
cat("Data files:\n")
for (file in c("hesim_setup.RData", "hesim_parameters.RData", 
               "hesim_transitions.RData", "hesim_costs_utilities.RData")) {
  if (file.exists(file.path("data", file))) {
    cat("  ✓ data/", file, "\n")
  } else {
    cat("  ✗ data/", file, " (missing)\n")
  }
}

cat("\nResult files:\n") 
if (file.exists("output/results/hesim_phase1_results.RData")) {
  cat("  ✓ output/results/hesim_phase1_results.RData\n")
} else {
  cat("  ✗ output/results/hesim_phase1_results.RData (missing)\n") 
}

cat("\nFigure files:\n")
for (file in c("hesim_death_rates.png", "hesim_cost_effectiveness.png")) {
  if (file.exists(file.path("output/figures", file))) {
    cat("  ✓ output/figures/", file, "\n")
  } else {
    cat("  ✗ output/figures/", file, " (missing)\n")
  }
}

cat("\n=== NEXT STEPS ===\n")
cat("1. Review results in output/results/ and output/figures/\n")
cat("2. If rates are too low, adjust risk_strata parameters in 01-model-setup.R\n")
cat("3. For Phase 2, expand to full 1000 risk strata and add PSA\n")
cat("4. For Phase 3, add risk prediction accuracy analysis\n")

cat("\n🎉 Hesim Phase 1 implementation complete! 🎉\n")