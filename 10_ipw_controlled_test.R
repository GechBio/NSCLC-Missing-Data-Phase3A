# ==============================================================================
# Script 10: IPW Controlled Test Under a Pure-MAR Design
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Isolate the source of IPW's bias in the main Monte Carlo. The main
#           simulation retains the original Unknown patients (fixed MNAR
#           component) and injects a second mechanism (MCAR / MAR / MNAR) on
#           top. Under a pure-MAR injection, IPW should recover nominal
#           performance IF the fixed MNAR component were absent. This script
#           removes the original Unknown patients entirely and runs a pure-MAR
#           injection, providing the positive control referenced in
#           Results §3.6 (Supplementary Table S3).
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Setup & Pure-MAR Cohort
#   MODULE 2: Controlled Monte Carlo
#   MODULE 3: Serialization & Export
#
# Save as:  10_ipw_controlled_test.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/Phase_3A/
#
# Input:    phase3a_base.rds
# Output:   TableS3_IPW_Controlled_Test.csv
#
# Depends:  simulation_helpers.R              (helper library)
#           05_simulation_monte_carlo.R       (sourced without execution)
#
# Seeds:    per-iteration seed = 123, n_iter = 100
# ==============================================================================


# ##############################################################################
# MODULE 1: Setup & Pure-MAR Cohort
# ##############################################################################

# --- Block 1.1: Load packages and construct pure-MAR cohort ---
init_ipw_control <- function() {
  cat("\n>>> MODULE 1, BLOCK 1.1: Loading Pure-MAR Cohort\n")
  
  if (!file.exists("phase3a_base.rds")) {
    stop("Run 00_base_setup_phase3a.R first.")
  }
  suppressPackageStartupMessages({
    library(survival)
    library(dplyr)
  })
  
  if (!file.exists("simulation_helpers.R")) {
    stop("Missing simulation_helpers.R (helper library for Script 10).")
  }
  source("simulation_helpers.R", local = FALSE)
  
  if (!file.exists("05_simulation_monte_carlo.R")) {
    stop("Missing 05_simulation_monte_carlo.R.")
  }
  old_opt <- getOption("phase3a.mc_execute")
  options(phase3a.mc_execute = FALSE)
  source("05_simulation_monte_carlo.R", local = FALSE)
  options(phase3a.mc_execute = old_opt)
  
  base   <- readRDS("phase3a_base.rds")
  cohort <- base$cohort_B_expr
  
  # Remove the original Unknown patients to eliminate the fixed MNAR component
  cohort <- cohort[cohort$stage_group != "Unknown", ]
  cohort$stage_group <- droplevels(factor(cohort$stage_group,
                                          levels = c("Early", "Advanced")))
  cohort
}


# ##############################################################################
# MODULE 2: Controlled Monte Carlo
# ##############################################################################

# --- Block 2.1: Run 100 pure-MAR iterations ---
run_pure_mar_monte_carlo <- function(cohort, n_iter = 100) {
  cat(">>> MODULE 2, BLOCK 2.1: Running Pure-MAR Monte Carlo\n")
  
  gt <- cohort |>
    assign_latent_stages_param(unk_hr = 3.22, seed = 42) |>
    simulate_survival_param(adv_hr = 1.86, early_median = 1984, seed = 999)
  
  true_log_hr <- log(1.86)
  
  sims <- list()
  set.seed(123)
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  for (i in seq_len(n_iter)) {
    setTxtProgressBar(pb, i)
    d <- gen_mar(gt)
    for (s in c("S1", "S2", "S3", "S4")) {
      r <- evaluate_strategy(d, s)
      sims[[length(sims) + 1]] <- data.frame(
        Iteration = i, Strategy = s,
        Beta_Adv  = r$beta_adv, SE_Adv = r$se_adv,
        stringsAsFactors = FALSE
      )
    }
  }
  close(pb)
  
  raw <- do.call(rbind, sims)
  
  raw |>
    group_by(Strategy) |>
    summarize(
      Mean_Beta = mean(Beta_Adv, na.rm = TRUE),
      Bias      = mean(Beta_Adv - true_log_hr, na.rm = TRUE),
      Coverage  = mean((Beta_Adv - 1.96 * SE_Adv <= true_log_hr) &
                         (Beta_Adv + 1.96 * SE_Adv >= true_log_hr), na.rm = TRUE),
      .groups = "drop"
    )
}


# ##############################################################################
# MODULE 3: Serialization & Export
# ##############################################################################

# --- Block 3.1: Print and save ---
export_ipw_control <- function(summ) {
  cat("\n>>> MODULE 3, BLOCK 3.1: Exporting IPW Controlled Test\n")
  
  cat("\n=== Pure-MAR controlled test (documented-stage-only cohort) ===\n")
  print(summ, row.names = FALSE)
  write.csv(summ, "TableS3_IPW_Controlled_Test.csv", row.names = FALSE)
  
  cat("\nInterpretation:\n")
  cat("  - If S4 (IPW) bias is near zero (|Bias| < 0.05), the two-component\n")
  cat("    interpretation is confirmed: IPW fails in the main simulation\n")
  cat("    because of the fixed MNAR component, not because of an IPW bug.\n")
  cat("  - If S4 bias remains large, investigate the IPW branch of\n")
  cat("    evaluate_strategy() before submission.\n")
  cat("Saved TableS3_IPW_Controlled_Test.csv\n")
}


# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_ipw_control_pipeline <- function() {
  cat("========================================================\n")
  cat("Starting Script 10: IPW Controlled Test (Pure MAR)\n")
  cat("========================================================\n")
  
  cohort <- init_ipw_control()
  summ   <- run_pure_mar_monte_carlo(cohort, n_iter = 100)
  export_ipw_control(summ)
  
  cat("\n========================================================\n")
  cat("Script 10 Complete!\n")
  cat("========================================================\n")
}

# Execute
run_ipw_control_pipeline()