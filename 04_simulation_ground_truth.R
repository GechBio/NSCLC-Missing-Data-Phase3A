# ==============================================================================
# Script 04: Phase 3A Simulation Arm 
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
#
# Purpose:  Construct the calibrated "Ground Truth" reference population for 
#           the simulation arm.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Setup & Data Loading
#   MODULE 2: Latent Ground Truth Generation
#   MODULE 3: Survival Simulation Engine
#   MODULE 4: Missingness Injection (Test MNAR)
#   MODULE 5: Evaluation & Export
#
# Input:    phase3a_base.rds
# Output:   phase3a_sim_ground_truth.rds
# ==============================================================================

# ##############################################################################
# MODULE 1: Setup & Data Loading
# ##############################################################################

init_sim_setup <- function() {
  cat("\n>>> MODULE 1: Loading Data and Packages\n")
  
  if (!file.exists("phase3a_base.rds")) stop("Run 00_base_setup_phase3a.R first!")
  base_data <- readRDS("phase3a_base.rds")
  
  library(survival)
  library(dplyr)
  
  # We use the full N=867 cohort (Strategy B)
  cohort <- base_data$cohort_B_expr
  
  # Ensure stage_group is a factor with 'Early' as the reference
  cohort$stage_group <- factor(cohort$stage_group, levels = c("Early", "Advanced", "Unknown"))
  
  return(list(base_data = base_data, cohort = cohort))
}


# ##############################################################################
# MODULE 2: Latent Ground Truth Generation
# ##############################################################################

assign_latent_stages <- function(cohort) {
  cat(">>> MODULE 2: Assigning Latent True Stages\n")
  
  set.seed(42)
  n_unknown <- sum(cohort$stage_group == "Unknown")
  
  # We assign a continuous latent log-hazard to the Unknown patients.
  # To achieve a composite HR of 3.22, the mean log-hazard is log(3.22) = 1.169.
  # We add a small SD (0.2) to represent heterogeneity (mixture of Advanced/Ultra-High)
  latent_log_hr <- rnorm(n_unknown, mean = log(3.22), sd = 0.2)
  
  cohort$latent_severity <- 0
  cohort$latent_severity[cohort$stage_group == "Unknown"] <- latent_log_hr
  
  # Create the true_stage column for simulation purposes
  cohort$true_stage <- as.character(cohort$stage_group)
  cohort$true_stage[cohort$stage_group == "Unknown"] <- "Ultra_High"
  cohort$true_stage <- factor(cohort$true_stage, levels = c("Early", "Advanced", "Ultra_High"))
  
  return(cohort)
}


# ##############################################################################
# MODULE 3: Survival Simulation Engine
# ##############################################################################

simulate_survival <- function(cohort) {
  cat(">>> MODULE 3: Simulating Survival Times\n")
  
  # Baseline Weibull hazard
  # We calibrate lambda so the Early group median OS is approx 1984 days
  target_early_median <- 1984
  shape <- 1.2
  lambda <- log(2) / (target_early_median^shape)
  
  # Linear predictor (log hazard ratio)
  lp <- rep(0, nrow(cohort))
  
  # True effect of Advanced stage (log(1.86) = 0.62)
  lp[cohort$true_stage == "Advanced"] <- log(1.86)
  
  # True effect of Unknown stage (latent continuous)
  lp[cohort$true_stage == "Ultra_High"] <- cohort$latent_severity[cohort$true_stage == "Ultra_High"]
  
  # Simulate survival times
  U <- runif(nrow(cohort))
  T_sim <- (-log(U) / (lambda * exp(lp)))^(1/shape)
  
  # Admin censoring: Uniform(0, 5000) to ensure high mortality in the Unknown group
  cens_time <- runif(nrow(cohort), min = 0, max = 5000) 
  
  cohort$sim_time  <- pmin(T_sim, cens_time)
  cohort$sim_event <- as.integer(T_sim <= cens_time)
  
  return(cohort)
}

calc_calibration_metrics <- function(sim_data) {
  # Fit Cox model using *observed* stage labels
  fit <- coxph(Surv(sim_time, sim_event) ~ stage_group, data = sim_data)
  hr_unknown <- exp(coef(fit)["stage_groupUnknown"])
  hr_adv     <- exp(coef(fit)["stage_groupAdvanced"])
  
  # Unknown Mortality & Median OS
  unknown_data <- sim_data[sim_data$stage_group == "Unknown", ]
  mortality_rate <- mean(unknown_data$sim_event)
  
  surv_fit <- survfit(Surv(sim_time, sim_event) ~ 1, data = unknown_data)
  median_os <- summary(surv_fit)$table["median"]
  
  return(list(
    HR_Unknown = hr_unknown,
    HR_Advanced = hr_adv,
    Mortality_Unknown = mortality_rate,
    Median_OS_Unknown = median_os
  ))
}


# ##############################################################################
# MODULE 4: Missingness Injection (Test MNAR)
# ##############################################################################

inject_mnar <- function(sim_data, target_rate = 0.142) {
  cat("\n>>> MODULE 4: Injecting Test MNAR Missingness\n")
  
  # MNAR mechanism: probability depends on survival time and latent severity
  gamma_1 <- -0.5  
  gamma_2 <- 1.0   
  
  # Find gamma_0 such that overall missingness is 14.2%
  find_gamma_0 <- function(g0) {
    logit_p <- g0 + gamma_1 * log(sim_data$sim_time + 1) + gamma_2 * sim_data$latent_severity
    p <- 1 / (1 + exp(-logit_p))
    mean(p) - target_rate
  }
  
  gamma_0 <- uniroot(find_gamma_0, c(-20, 20))$root
  
  logit_p <- gamma_0 + gamma_1 * log(sim_data$sim_time + 1) + gamma_2 * sim_data$latent_severity
  p_missing <- 1 / (1 + exp(-logit_p))
  
  sim_data$observed_stage <- as.character(sim_data$stage_group)
  sim_data$observed_stage[runif(nrow(sim_data)) < p_missing] <- "Unknown"
  sim_data$observed_stage <- factor(sim_data$observed_stage, levels = c("Early", "Advanced", "Unknown"))
  
  cat("Injected missingness rate:", round(mean(sim_data$observed_stage == "Unknown"), 3), "\n")
  
  # Check observed composite HR after injection
  fit_obs <- coxph(Surv(sim_time, sim_event) ~ observed_stage, data = sim_data)
  hr_obs <- exp(coef(fit_obs)["observed_stageUnknown"])
  
  cat("Observed Composite Unknown HR after MNAR injection:", round(hr_obs, 2), "\n")
  
  return(sim_data)
}


# ##############################################################################
# MODULE 5: Evaluation & Export
# ##############################################################################

export_ground_truth <- function(cohort) {
  cat("\n>>> MODULE 5: Exporting Ground Truth\n")
  
  set.seed(999)
  gt_data <- simulate_survival(cohort)
  
  metrics <- calc_calibration_metrics(gt_data)
  cat("\n--- Final Ground Truth Metrics ---\n")
  cat("Advanced HR (Target 1.86):", round(metrics$HR_Advanced, 2), "\n")
  cat("Unknown Composite HR (Target 3.22):", round(metrics$HR_Unknown, 2), "\n")
  cat("Unknown Mortality (Target 0.943):", round(metrics$Mortality_Unknown, 3), "\n")
  cat("Unknown Median OS (Target 692):", round(metrics$Median_OS_Unknown, 0), "days\n")
  
  # Test MNAR Injection on the SIMULATED data (not the latent cohort)
  test_mnar_data <- inject_mnar(gt_data)
  
  saveRDS(list(
    ground_truth_data = gt_data,
    final_metrics = metrics
  ), "phase3a_sim_ground_truth.rds")
  
  cat("\nSaved phase3a_sim_ground_truth.rds successfully.\n")
}


# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_simulation_ground_truth <- function() {
  cat("========================================================\n")
  cat("Starting Phase 3A Simulation Ground Truth Calibration\n")
  cat("========================================================\n")
  
  setup <- init_sim_setup()
  cohort_latent <- assign_latent_stages(setup$cohort)
  export_ground_truth(cohort_latent)
  
  cat("\n========================================================\n")
  cat("Ground Truth Calibration Complete!\n")
  cat("========================================================\n")
}

# Execute
run_simulation_ground_truth()