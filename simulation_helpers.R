# ==============================================================================
# Script 04R: Phase 3A Simulation Arm - Parameterized Ground Truth Helper
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Parameterized versions of the latent-stage assignment and survival
#           simulation functions from Script 04. Defaults reproduce the original
#           manuscript targets exactly. The parameterization exists so that
#           Script 10 (calibration sensitivity) and Script 11 (IPW controlled
#           test) can vary the four calibration targets without duplicating the
#           simulation engine.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Parameterized Latent-Stage Assignment
#   MODULE 2: Parameterized Survival Simulation
#   MODULE 3: Achieved-Metrics Reporter
#
# Save as:  simulation_helpers.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/Phase_3A/
# Role:     Helper library, not a pipeline step. Sourced by Scripts 10 and 11.
#           The original Script 04 (Step 4) remains the pipeline step that
#           produces phase3a_sim_ground_truth.rds.
#
# Input:    none (functions only)
# Output:   none (functions only; does NOT write files)
#
# NOTE:     This file is a helper, not a pipeline step. It is sourced by
#           Script 10 and Script 11. It is never sourced by 08_master_pipeline.R.
# ==============================================================================


# ##############################################################################
# MODULE 1: Parameterized Latent-Stage Assignment
# ##############################################################################

# --- Block 1.1: assign_latent_stages_param() ---
assign_latent_stages_param <- function(cohort,
                                       unk_hr       = 3.22,
                                       unk_sd_loghr = 0.2,
                                       seed         = 42) {
  set.seed(seed)
  n_unknown <- sum(cohort$stage_group == "Unknown")
  
  latent_log_hr <- rnorm(n_unknown, mean = log(unk_hr), sd = unk_sd_loghr)
  
  cohort$latent_severity <- 0
  cohort$latent_severity[cohort$stage_group == "Unknown"] <- latent_log_hr
  
  cohort$true_stage <- as.character(cohort$stage_group)
  cohort$true_stage[cohort$stage_group == "Unknown"] <- "Ultra_High"
  cohort$true_stage <- factor(cohort$true_stage,
                              levels = c("Early", "Advanced", "Ultra_High"))
  
  cohort
}


# ##############################################################################
# MODULE 2: Parameterized Survival Simulation
# ##############################################################################

# --- Block 2.1: simulate_survival_param() ---
simulate_survival_param <- function(cohort,
                                    adv_hr       = 1.86,
                                    early_median = 1984,
                                    shape        = 1.2,
                                    censor_max   = 5000,
                                    seed         = 999) {
  set.seed(seed)
  
  lambda <- log(2) / (early_median^shape)
  lp     <- rep(0, nrow(cohort))
  
  lp[cohort$true_stage == "Advanced"]   <- log(adv_hr)
  lp[cohort$true_stage == "Ultra_High"] <-
    cohort$latent_severity[cohort$true_stage == "Ultra_High"]
  
  U       <- runif(nrow(cohort))
  T_sim   <- (-log(U) / (lambda * exp(lp)))^(1 / shape)
  cens    <- runif(nrow(cohort), min = 0, max = censor_max)
  
  cohort$sim_time  <- pmin(T_sim, cens)
  cohort$sim_event <- as.integer(T_sim <= cens)
  
  cohort
}


# ##############################################################################
# MODULE 3: Achieved-Metrics Reporter
# ##############################################################################

# --- Block 3.1: report_calibration() ---
report_calibration <- function(sim_data) {
  fit <- coxph(Surv(sim_time, sim_event) ~ stage_group, data = sim_data)
  
  unk <- sim_data[sim_data$stage_group == "Unknown", ]
  km  <- survfit(Surv(sim_time, sim_event) ~ 1, data = unk)
  
  data.frame(
    Advanced_HR    = round(exp(coef(fit)["stage_groupAdvanced"]), 3),
    Unknown_HR     = round(exp(coef(fit)["stage_groupUnknown"]),  3),
    Unknown_Mort   = round(mean(unk$sim_event), 4),
    Unknown_KM_Med = round(summary(km)$table["median"], 1),
    stringsAsFactors = FALSE
  )
}