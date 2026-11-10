# ==============================================================================
# Script 08: Phase 3A Master Execution Pipeline
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Execute the entire Phase 3A pipeline sequentially from start to 
#           finish. This script guarantees 100% reproducibility of all 
#           manuscript results, tables, and figures for the Zenodo repository.
#
# ARCHITECTURE: SEQUENTIAL EXECUTION
#   - Step 1: Base Setup & Data Assembly (Script 00)
#   - Step 2: Empirical Arm Modeling (Scripts 01, 02)
#   - Step 3: Empirical Figures (Scripts 03, 03b)
#   - Step 4: Simulation Arm (Scripts 04, 04b, 05, 06)
#   - Step 5: Reporting Checklist (Script 07)
#
# Save as:  08_master_pipeline.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/Phase_3A/
#
# Input:    All raw data (from parent directory) and prior scripts (00-07)
# Output:   All intermediate .rds files, Figures 2-6, and Table 5 (CSV)
#
# Dependencies: Requires scripts 00 through 07 to be in the same directory
# ==============================================================================

cat("========================================================\n")
cat("Starting Full Phase 3A Reproduction Pipeline\n")
cat("========================================================\n")

# --- Step 1: Base Setup & Data Assembly (N=867 Expr, N=701 Meth) ---
cat("\n>>> EXECUTING: Script 00 (Base Setup)\n")
source("00_base_setup_phase3a.R")

# --- Step 2: Empirical Arm Modeling ---
cat("\n>>> EXECUTING: Script 01 (Expression Modeling)\n")
source("01_empirical_expression.R")

cat("\n>>> EXECUTING: Script 02 (Methylation Modeling)\n")
source("02_empirical_methylation.R")

# --- Step 3: Empirical Figures ---
cat("\n>>> EXECUTING: Script 03 (C-index & HR Figures)\n")
source("03_generate_empirical_figures.R")

cat("\n>>> EXECUTING: Script 03b (Feature Selection Dot Plot)\n")
source("03b_feature_selection_upset.R")

# --- Step 4: Simulation Arm ---
cat("\n>>> EXECUTING: Script 04 (Simulation Ground Truth)\n")
source("04_simulation_ground_truth.R")

cat("\n>>> EXECUTING: Script 04b (Ground Truth KM Figure)\n")
source("04b_generate_simulation_ground_truth_figure.R")

cat("\n>>> EXECUTING: Script 05 (Monte Carlo Engine)\n")
source("05_simulation_monte_carlo.R")

cat("\n>>> EXECUTING: Script 06 (Combined Monte Carlo Figure)\n")
source("06_simulation_figures.R")

# --- Step 5: Reporting Checklist ---
cat("\n>>> EXECUTING: Script 07 (Reporting Checklist)\n")
source("07_reporting_checklist.R")

cat("\n========================================================\n")
cat("Full Phase 3A Reproduction Pipeline Complete!\n")
cat("All tables and figures have been generated successfully.\n")
cat("========================================================\n")