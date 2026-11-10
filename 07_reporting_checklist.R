# ==============================================================================
# Script 07: Phase 3A - Generate Reporting Checklist
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Output the final reporting checklist for prognostic biomarker 
#           studies with missing clinical covariates. This fulfills the 
#           manuscript's promise of providing a practical standardization guide.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Checklist Definition
#       - Block 1.1: Define the 7-point reporting framework
#   MODULE 2: Serialization & Export
#       - Block 2.1: Print to console and save as CSV
#
# Save as:  07_reporting_checklist.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/Phase_3A/
#
# Input:    None (Standalone definition)
# Output:   Table_Reporting_Checklist.csv
#
# Dependencies: base R (no external packages required)
# ==============================================================================

# ##############################################################################
# MODULE 1: Checklist Definition
# ##############################################################################

# --- Block 1.1: Define the 7-point reporting framework ---
define_checklist <- function() {
  cat("\n>>> MODULE 1, BLOCK 1.1: Defining Checklist\n")
  
  checklist <- data.frame(
    Item = 1:7,
    Domain = c("Reporting", "Reporting", "Analysis", "Analysis", "Sensitivity", "Sensitivity", "Interpretation"),
    Recommendation = c(
      "Report the proportion of missingness for each clinical covariate.",
      "Report the distribution of missingness by key demographic and histologic strata.",
      "Test whether missingness is associated with the outcome; if so, treat it as informative.",
      "Pre-specify the primary missingness-handling strategy and justify it based on the likely mechanism.",
      "Report sensitivity analyses across at least two alternative handling strategies (e.g., complete-case and missing-retained).",
      "Report feature-selection stability under the primary strategy and at least one alternative.",
      "Avoid claiming unbiased discovery when complete-case analysis has been applied to a cohort with informative missingness."
    )
  )
  
  return(checklist)
}

# ##############################################################################
# MODULE 2: Serialization & Export
# ##############################################################################

# --- Block 2.1: Print to console and save as CSV ---
export_checklist <- function(checklist) {
  cat(">>> MODULE 2, BLOCK 2.1: Exporting Checklist\n")
  
  cat("\n========================================================\n")
  cat("Phase 3A Reporting Checklist\n")
  cat("========================================================\n")
  print(checklist)
  
  write.csv(checklist, "Table_Reporting_Checklist.csv", row.names = FALSE)
  cat("\nSaved Table_Reporting_Checklist.csv successfully.\n")
}

# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_checklist_pipeline <- function() {
  cat("========================================================\n")
  cat("Starting Phase 3A Checklist Generation\n")
  cat("========================================================\n")
  
  my_checklist <- define_checklist()
  export_checklist(my_checklist)
  
  cat("\n========================================================\n")
  cat("Checklist Generation Complete!\n")
  cat("========================================================\n")
}

# Execute
run_checklist_pipeline()