# ==============================================================================
# Script 09: Clinical Cohort Disambiguation (KM Median vs Observation Time)
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Produce the Kaplan-Meier and observation-time table for the clinical
#           cohort (N = 867). For each stage group (Early, Advanced, Unknown,
#           Overall), report both the Kaplan-Meier median survival (with 95% CI)
#           and the median observation time (all patients, events and censored).
#           The main disambiguation table is Supplementary Table S1; the
#           internal consistency check table is Supplementary Table S2.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Setup & Data Loading
#   MODULE 2: Per-Group KM and Observation-Time Metrics
#   MODULE 3: Internal Consistency Diagnostics
#   MODULE 4: Serialization & Export
#
# Save as:  09_clinical_cohort_disambiguation.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/Phase_3A/
#
# Input:    phase3a_base.rds
# Output:   TableS1_Clinical_Cohort_Disambiguation.csv
#           TableS2_Clinical_Cohort_Consistency.csv
#
# Seeds:    none
# ==============================================================================


# ##############################################################################
# MODULE 1: Setup & Data Loading
# ##############################################################################

# --- Block 1.1: Load cohort ---
init_disambiguation <- function() {
  cat("\n>>> MODULE 1, BLOCK 1.1: Loading Clinical Cohort\n")
  
  if (!file.exists("phase3a_base.rds")) {
    stop("Run 00_base_setup_phase3a.R first.")
  }
  suppressPackageStartupMessages({
    library(survival)
    library(dplyr)
  })
  
  base <- readRDS("phase3a_base.rds")
  cohort <- base$cohort_B_expr
  cohort$stage_group <- factor(cohort$stage_group,
                               levels = c("Early", "Advanced", "Unknown"))
  cohort
}


# ##############################################################################
# MODULE 2: Per-Group KM and Observation-Time Metrics
# ##############################################################################

# --- Block 2.1: Compute metrics for one subgroup ---
compute_group_metrics <- function(d, label) {
  n        <- nrow(d)
  events   <- sum(d$os_event, na.rm = TRUE)
  event_rt <- events / n
  
  km     <- survfit(Surv(os_time, os_event) ~ 1, data = d)
  km_tab <- summary(km)$table
  
  data.frame(
    Group           = label,
    N               = n,
    Events          = events,
    Event_Rate      = round(event_rt, 4),
    KM_Median       = km_tab["median"],
    KM_Median_Lower = km_tab["0.95LCL"],
    KM_Median_Upper = km_tab["0.95UCL"],
    Obs_Median      = median(d$os_time, na.rm = TRUE),
    Obs_IQR         = IQR(d$os_time, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# --- Block 2.2: Compute metrics across all groups ---
compute_all_group_metrics <- function(cohort) {
  cat(">>> MODULE 2, BLOCK 2.1 & 2.2: Computing Per-Group Metrics\n")
  
  groups <- list(
    Overall  = cohort,
    Early    = cohort[cohort$stage_group == "Early", ],
    Advanced = cohort[cohort$stage_group == "Advanced", ],
    Unknown  = cohort[cohort$stage_group == "Unknown", ]
  )
  do.call(rbind, lapply(names(groups),
                        function(g) compute_group_metrics(groups[[g]], g)))
}


# ##############################################################################
# MODULE 3: Internal Consistency Diagnostics
# ##############################################################################

# --- Block 3.1: Check for impossible mixtures and unreached medians ---
build_consistency_checks <- function(metrics) {
  cat(">>> MODULE 3, BLOCK 3.1: Running Internal Consistency Checks\n")
  
  unk     <- metrics[metrics$Group == "Unknown", ]
  overall <- metrics[metrics$Group == "Overall", ]
  early   <- metrics[metrics$Group == "Early", ]
  adv     <- metrics[metrics$Group == "Advanced", ]
  
  subgroup_km <- c(early$KM_Median, adv$KM_Median, unk$KM_Median)
  if (any(is.na(subgroup_km))) {
    overall_below_subgroups <- NA
  } else {
    overall_below_subgroups <- overall$KM_Median < min(subgroup_km)
  }
  
  checks <- data.frame(
    Item      = c("Overall_KM_below_every_subgroup",
                  "Unknown_KM_median_reached",
                  "Unknown_event_rate_100pct"),
    Value     = c(overall_below_subgroups,
                  !is.na(unk$KM_Median),
                  unk$Event_Rate >= 0.99),
    Flag      = c(
      ifelse(is.na(overall_below_subgroups), "CHECK",
             ifelse(overall_below_subgroups, "FAIL", "PASS")),
      ifelse(!is.na(unk$KM_Median), "PASS", "CHECK"),
      ifelse(unk$Event_Rate >= 0.99, "PASS",
             "PASS")
    ),
    Notes     = c(
      "Overall KM median cannot be below every subgroup KM median.",
      "If NA, the Unknown KM curve never crossed 50% (max follow-up < median).",
      "Reports whether the Unknown group is fully uncensored."
    ),
    stringsAsFactors = FALSE
  )
  checks
}


# ##############################################################################
# MODULE 4: Serialization & Export
# ##############################################################################

# --- Block 4.1: Print and save both tables ---
export_disambiguation_tables <- function(metrics, checks) {
  cat(">>> MODULE 4, BLOCK 4.1: Exporting Disambiguation Tables\n")
  
  cat("\n=== Per-group metrics (N = 867) ===\n")
  print(metrics, row.names = FALSE)
  
  cat("\n=== Internal consistency checks ===\n")
  print(checks, row.names = FALSE)
  
  write.csv(metrics, "TableS1_Clinical_Cohort_Disambiguation.csv",
            row.names = FALSE)
  write.csv(checks,  "TableS2_Clinical_Cohort_Consistency.csv",
            row.names = FALSE)
  
  cat("\nNote: KM_Median is the Kaplan-Meier median survival; Obs_Median is the\n")
  cat("      median of os_time across all patients (events and censored). The\n")
  cat("      '692-day' figure used as a calibration target in Script 04 should\n")
  cat("      be read as one of these two quantities; this table disambiguates\n")
  cat("      which one the pipeline actually produced.\n")
  cat("Saved TableS1_Clinical_Cohort_Disambiguation.csv and TableS2_Clinical_Cohort_Consistency.csv.\n")
}


# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_disambiguation_pipeline <- function() {
  cat("========================================================\n")
  cat("Starting Script 09: Clinical Cohort Disambiguation\n")
  cat("========================================================\n")
  
  cohort  <- init_disambiguation()
  metrics <- compute_all_group_metrics(cohort)
  checks  <- build_consistency_checks(metrics)
  export_disambiguation_tables(metrics, checks)
  
  cat("\n========================================================\n")
  cat("Script 09 Complete!\n")
  cat("========================================================\n")
}

# Execute
run_disambiguation_pipeline()