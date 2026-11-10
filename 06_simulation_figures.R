# ==============================================================================
# Script 06: Phase 3A Simulation Arm - Combined Monte Carlo Figure
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Generate a single, comprehensive Dot-and-Whisker plot that replaces 
#           the separate Bias and Coverage bar charts. This shows the distribution 
#           of the 100 Monte Carlo estimates, the mean estimate, the 95% CI of the 
#           mean, and visually indicates whether each iteration covered the true 
#           hazard ratio.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load Simulation Results
#   MODULE 2: Data Processing
#       - Block 2.1: Calculate Coverage and Summary Statistics
#   MODULE 3: Plot Generation & Export
#       - Block 3.1: Generate and Save Combined Plot
#
# Save as:  06_simulation_figures_combined.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/Phase_3A/
#
# Input:    phase3a_sim_monte_carlo_results.rds
# Output:   Figure6_Simulation_Combined.pdf
#
# Dependencies: ggplot2, dplyr
# ==============================================================================

# ##############################################################################
# MODULE 1: Environment & Dependency Management
# ##############################################################################

# --- Block 1.1: Load Required Packages ---
init_combined_setup <- function() {
  cat("\n>>> MODULE 1, BLOCK 1.1: Loading Packages\n")
  
  library(ggplot2)
  library(dplyr)
}

# --- Block 1.2: Load Simulation Results ---
load_simulation_data <- function() {
  cat(">>> MODULE 1, BLOCK 1.2: Loading Simulation Results\n")
  
  if (!file.exists("phase3a_sim_monte_carlo_results.rds")) {
    stop("Run Script 05 first to generate the simulation results.")
  }
  
  sim_results <- readRDS("phase3a_sim_monte_carlo_results.rds")
  raw_df <- sim_results$raw_results
  
  # True log hazard ratio (log(1.86))
  true_log_hr <- log(1.86)
  
  return(list(raw_df = raw_df, true_log_hr = true_log_hr))
}

# ##############################################################################
# MODULE 2: Data Processing
# ##############################################################################

# --- Block 2.1: Calculate Coverage and Summary Statistics ---
process_plot_data <- function(raw_df, true_log_hr) {
  cat(">>> MODULE 2, BLOCK 2.1: Calculating Coverage and Summary Statistics\n")
  
  # Calculate if the 95% CI of each iteration covered the true value
  raw_df <- raw_df %>%
    mutate(
      lower = Beta_Adv - 1.96 * SE_Adv,
      upper = Beta_Adv + 1.96 * SE_Adv,
      covered = lower <= true_log_hr & upper >= true_log_hr,
      Strategy = factor(Strategy, levels = c("S1", "S2", "S3", "S4"),
                        labels = c("S1: Complete-Case", "S2: Unknown Retention",
                                   "S3: Multiple Imputation", "S4: IPW"))
    )
  
  # Calculate summary statistics (mean and 95% CI of the mean)
  summary_df <- raw_df %>%
    group_by(Mechanism, Strategy) %>%
    summarize(
      mean_beta = mean(Beta_Adv, na.rm = TRUE),
      sd_beta = sd(Beta_Adv, na.rm = TRUE),
      n = n(),
      se_beta = sd_beta / sqrt(n),
      ci_lower = mean_beta - 1.96 * se_beta,
      ci_upper = mean_beta + 1.96 * se_beta,
      .groups = "drop"
    )
  
  return(list(raw_df = raw_df, summary_df = summary_df))
}

# ##############################################################################
# MODULE 3: Plot Generation & Export
# ##############################################################################

# --- Block 3.1: Generate and Save Combined Plot ---
generate_combined_plot <- function(processed_data, true_log_hr) {
  cat(">>> MODULE 3, BLOCK 3.1: Generating Combined Monte Carlo Plot\n")
  
  raw_df <- processed_data$raw_df
  summary_df <- processed_data$summary_df
  
  p <- ggplot(raw_df, aes(x = Strategy, y = Beta_Adv)) +
    # True value line
    geom_hline(yintercept = true_log_hr, linetype = "dashed", color = "red", linewidth = 1) +
    
    # Jittered raw points (colored by whether they covered the truth)
    geom_jitter(aes(color = covered), width = 0.2, alpha = 0.4, size = 1.5) +
    
    # Summary point range (Mean and 95% CI of the mean) - Removed 'fatten'
    geom_pointrange(data = summary_df, 
                    aes(y = mean_beta, ymin = ci_lower, ymax = ci_upper),
                    color = "black", size = 0.8, linewidth = 0.8) +
    
    facet_wrap(~ Mechanism) +
    
    # Colors for coverage (Blue = Covered, Red = Not Covered)
    scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"),
                       labels = c("TRUE" = "95% CI Covered Truth", "FALSE" = "95% CI Missed Truth")) +
    
    labs(
      title = "Monte Carlo Simulation: Estimated vs. True Advanced Stage Log HR",
      subtitle = paste0("True Log HR = ", round(true_log_hr, 2), " (Red Dashed Line). Black points show mean and 95% CI of the estimates."),
      x = "Handling Strategy",
      y = "Estimated Log Hazard Ratio",
      color = "Iteration Coverage"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 15, hjust = 1)
    )
  
  ggsave("Figure6_Simulation_Combined.pdf", p, width = 10, height = 6)
  cat("Saved Figure6_Simulation_Combined.pdf successfully.\n")
}

# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_combined_pipeline <- function() {
  cat("========================================================\n")
  cat("Starting Combined Monte Carlo Figure Generation\n")
  cat("========================================================\n")
  
  init_combined_setup()
  setup <- load_simulation_data()
  processed_data <- process_plot_data(setup$raw_df, setup$true_log_hr)
  generate_combined_plot(processed_data, setup$true_log_hr)
  
  cat("\n========================================================\n")
  cat("Combined Figure Generation Complete!\n")
  cat("========================================================\n")
}

# Execute
run_combined_pipeline()