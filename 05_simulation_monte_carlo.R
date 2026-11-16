# ==============================================================================
# Script 05: Phase 3A Simulation Arm - Monte Carlo Engine 
#
# Author:   Getachew Derbew Belay
# Date:     November 2026 (revised)
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Run the Monte Carlo simulation to benchmark missing-data handling 
#           strategies (S1-S4) under MCAR, MAR, and MNAR mechanisms. 
#           Calculates bias, coverage, and MSE against the known ground truth.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Load Required Packages
#       - Block 1.2: Load Ground Truth Data
#   MODULE 2: Missingness Generators (Two-Component Design)
#       - Block 2.1: MCAR Generator
#       - Block 2.2: MAR Generator
#       - Block 2.3: MNAR Generator
#   MODULE 3: Handling Strategies (S1, S2, S3, S4)
#       - Block 3.1: Core Strategy Evaluation Function
#   MODULE 4: Monte Carlo Engine
#       - Block 4.1: Execution Loop (with Seed)
#   MODULE 5: Summarization & Export
#       - Block 5.1: Calculate Bias, MSE, Coverage
#       - Block 5.2: Save Results
#
# Save as:  05_simulation_monte_carlo.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/Phase_3A/
#
# Input:    phase3a_sim_ground_truth.rds
# Output:   phase3a_sim_monte_carlo_results.rds
#
# Dependencies: survival, dplyr
#
# REVISION NOTES (Nov 2026):
#   - S4 (IPW) weights corrected: was 1 / P(missing|X), now 1 / P(observed|X).
#     The previous direction systematically upweighted the "easy-to-observe"
#     subjects and downweighted the "hard-to-observe" ones, inflating bias in
#     all mechanisms. Fixed by computing weights as 1 / (1 - P(missing|X)).
#   - S4 weights are now stabilized (mean weight ~= 1) and truncated at the
#     99th percentile, matching the description in the manuscript ("stabilized,
#     truncated weights").
#   - S4 uses robust (sandwich) standard errors via cluster = row_id; the
#     model-based SE from weighted coxph underestimates variance for IPW.
#   - No change to S1, S2, S3. The two-component missingness design (fixed
#     MNAR + injected mechanism) is documented in Module 2.
# ==============================================================================

# ##############################################################################
# MODULE 1: Environment & Dependency Management
# ##############################################################################

# --- Block 1.1: Load Required Packages ---
init_environment <- function() {
  cat("\n>>> MODULE 1, BLOCK 1.1: Loading Packages\n")
  
  packages <- c("survival", "dplyr")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# --- Block 1.2: Load Ground Truth Data ---
load_ground_truth <- function() {
  cat(">>> MODULE 1, BLOCK 1.2: Loading Ground Truth Data\n")
  
  if (!file.exists("phase3a_sim_ground_truth.rds")) {
    stop("Run Script 04 first to generate the ground truth.")
  }
  
  gt_data <- readRDS("phase3a_sim_ground_truth.rds")
  cohort <- gt_data$ground_truth_data
  
  true_beta_adv <- log(1.86)
  true_beta_unknown <- log(3.22) 
  
  return(list(cohort = cohort, 
              true_beta_adv = true_beta_adv, 
              true_beta_unknown = true_beta_unknown))
}


# ##############################################################################
# MODULE 2: Missingness Generators (Two-Component Design)
# ##############################################################################

# -----------------------------------------------------------------------------
# TWO-COMPONENT MISSINGNESS DESIGN (see Methods §2.6 of manuscript)
# -----------------------------------------------------------------------------
# The `cohort` passed to gen_*() already contains the ORIGINAL Unknown patients
# (n = 123, ~14.2% of N = 867). Their missingness is driven by unobserved
# latent severity, making it MNAR by construction. On top of this FIXED MNAR
# component, the functions below inject an ADDITIONAL 14.2% marginal
# probability of missingness under MCAR, MAR, or MNAR. The `observed_stage`
# column therefore mixes both sources:
#   (1) original Unknown (fixed MNAR) — always present
#   (2) injected Unknown (mechanism-specific) — added on top
# This mirrors real-world TCGA data, where undocumented stage reflects both a
# genuinely high-risk phenotype and administrative incompleteness.
# Script 11 (verify_04 in earlier drafts) isolates the injected component by
# removing the original Unknown patients entirely, providing the positive
# control for S4.
# -----------------------------------------------------------------------------

# --- Block 2.1: MCAR Generator ---
gen_mcar <- function(data, target_rate = 0.142) {
  n <- nrow(data)
  drop_idx <- sample(1:n, size = round(n * target_rate))
  
  data$observed_stage <- as.character(data$stage_group)
  data$observed_stage[drop_idx] <- "Unknown"
  data$observed_stage <- factor(data$observed_stage, levels = c("Early", "Advanced", "Unknown"))
  return(data)
}

# --- Block 2.2: MAR Generator (Depends on Age & Histology) ---
gen_mar <- function(data, target_rate = 0.142) {
  find_intercept <- function(b0) {
    p <- 1 / (1 + exp(-(b0 + 0.02 * data$age + ifelse(data$histology == "LUSC", 0.5, 0))))
    mean(p) - target_rate
  }
  b0 <- uniroot(find_intercept, c(-20, 20))$root
  
  logit_p <- b0 + 0.02 * data$age + ifelse(data$histology == "LUSC", 0.5, 0)
  p_missing <- 1 / (1 + exp(-logit_p))
  
  data$observed_stage <- as.character(data$stage_group)
  data$observed_stage[runif(nrow(data)) < p_missing] <- "Unknown"
  data$observed_stage <- factor(data$observed_stage, levels = c("Early", "Advanced", "Unknown"))
  return(data)
}

# --- Block 2.3: MNAR Generator (Depends on Survival & Latent Severity) ---
gen_mnar <- function(data, target_rate = 0.142) {
  gamma_1 <- -0.5  
  gamma_2 <- 1.0   
  
  find_gamma_0 <- function(g0) {
    logit_p <- g0 + gamma_1 * log(data$sim_time + 1) + gamma_2 * data$latent_severity
    p <- 1 / (1 + exp(-logit_p))
    mean(p) - target_rate
  }
  gamma_0 <- uniroot(find_gamma_0, c(-20, 20))$root
  
  logit_p <- gamma_0 + gamma_1 * log(data$sim_time + 1) + gamma_2 * data$latent_severity
  p_missing <- 1 / (1 + exp(-logit_p))
  
  data$observed_stage <- as.character(data$stage_group)
  data$observed_stage[runif(nrow(data)) < p_missing] <- "Unknown"
  data$observed_stage <- factor(data$observed_stage, levels = c("Early", "Advanced", "Unknown"))
  return(data)
}


# ##############################################################################
# MODULE 3: Handling Strategies (S1, S2, S3, S4)
# ##############################################################################

# --- Block 3.1: Core Strategy Evaluation Function ---
evaluate_strategy <- function(data, strategy) {
  tryCatch({
    # S1: Complete-Case Deletion
    if (strategy == "S1") {
      dat_sub <- data[data$observed_stage != "Unknown", ]
      dat_sub$observed_stage <- droplevels(dat_sub$observed_stage)
      fit <- coxph(Surv(sim_time, sim_event) ~ observed_stage, data = dat_sub)
      beta_unknown <- NA 
      beta_adv <- coef(fit)["observed_stageAdvanced"]
      se_adv <- sqrt(vcov(fit)["observed_stageAdvanced", "observed_stageAdvanced"])
      
      # S2: Unknown Retention
    } else if (strategy == "S2") {
      dat_sub <- data
      fit <- coxph(Surv(sim_time, sim_event) ~ observed_stage, data = dat_sub)
      beta_unknown <- coef(fit)["observed_stageUnknown"]
      beta_adv <- coef(fit)["observed_stageAdvanced"]
      se_adv <- sqrt(vcov(fit)["observed_stageAdvanced", "observed_stageAdvanced"])
      
      # S3: Multiple Imputation (Manual, robust to crashes)
    } else if (strategy == "S3") {
      dat_sub <- data[, c("sim_time", "sim_event", "age", "histology", "observed_stage")]
      
      # Isolate observed and missing data
      dat_obs <- dat_sub[dat_sub$observed_stage != "Unknown", ]
      dat_miss <- dat_sub[dat_sub$observed_stage == "Unknown", ]
      
      if (nrow(dat_miss) > 0) {
        # 1. Fit imputation model on observed data (binary outcome: Advanced vs Early)
        model_data <- dat_obs
        model_data$is_advanced <- as.integer(model_data$observed_stage == "Advanced")
        imp_model <- glm(is_advanced ~ age + histology, data = model_data, family = binomial)
        
        # 2. Predict probabilities for missing patients
        pred_probs <- predict(imp_model, newdata = dat_miss, type = "response")
        
        # 3. Impute 20 times
        beta_adv_imps <- numeric(20)
        se_adv_imps <- numeric(20)
        
        for (m in 1:20) {
          # Generate imputed stages
          imputed_stage <- ifelse(runif(nrow(dat_miss)) < pred_probs, "Advanced", "Early")
          
          # Build completed dataset (clean rbind without extra columns)
          dat_comp <- rbind(dat_obs, dat_miss)
          dat_comp$observed_stage[dat_comp$observed_stage == "Unknown"] <- imputed_stage
          dat_comp$observed_stage <- factor(dat_comp$observed_stage, levels = c("Early", "Advanced"))
          
          # Fit Cox model
          fit <- coxph(Surv(sim_time, sim_event) ~ observed_stage, data = dat_comp)
          beta_adv_imps[m] <- coef(fit)["observed_stageAdvanced"]
          se_adv_imps[m] <- sqrt(vcov(fit)["observed_stageAdvanced", "observed_stageAdvanced"])
        }
        
        # 4. Pool using Rubin's Rules
        beta_adv <- mean(beta_adv_imps)
        se_adv <- sqrt(mean(se_adv_imps^2) + (1 + 1/20) * var(beta_adv_imps))
        
      } else {
        beta_adv <- NA
        se_adv <- NA
      }
      beta_unknown <- NA 
      
      # S4: Inverse Probability Weighting (CORRECTED — see revision notes)
    } else if (strategy == "S4") {
      dat_full <- data
      dat_full$missing <- as.integer(dat_full$observed_stage == "Unknown")
      
      # Propensity model: P(missing | age, histology)
      ps_model <- glm(missing ~ age + histology, data = dat_full, family = binomial)
      
      # Complete-case subset (drop Unknown patients)
      dat_sub <- dat_full[dat_full$missing == 0, ]
      dat_sub$observed_stage <- droplevels(dat_sub$observed_stage)
      
      # P(missing | X) for observed subjects; truncate propensities to [0.01, 0.99]
      p_missing <- predict(ps_model, newdata = dat_sub, type = "response")
      p_missing <- pmax(pmin(p_missing, 0.99), 0.01)
      
      # CORRECTED IPW: weight = 1 / P(observed | X) = 1 / (1 - P(missing | X))
      # (Previous version used 1 / P(missing | X), the inverse of what is needed.
      #  That systematically upweighted the easy-to-observe subjects.)
      # Stabilize by the marginal probability of being observed so weights
      # have mean ~= 1, reducing variance without changing the point estimate.
      p_obs_marginal <- mean(dat_full$missing == 0)
      dat_sub$weights <- p_obs_marginal / (1 - p_missing)
      
      # Truncate extreme weights at the 99th percentile (stability)
      w_cap <- as.numeric(quantile(dat_sub$weights, 0.99, na.rm = TRUE))
      dat_sub$weights <- pmin(dat_sub$weights, w_cap)
      
      # Weighted Cox with ROBUST (sandwich) standard errors.
      # The model-based SE from weighted coxph underestimates variance for IPW,
      # so we cluster on unique row IDs to obtain the robust variance.
      dat_sub$row_id <- seq_len(nrow(dat_sub))
      fit <- coxph(Surv(sim_time, sim_event) ~ observed_stage,
                   data = dat_sub, weights = weights, cluster = row_id)
      
      beta_unknown <- NA 
      beta_adv <- coef(fit)["observed_stageAdvanced"]
      se_adv <- sqrt(vcov(fit)["observed_stageAdvanced", "observed_stageAdvanced"])
    }
    
    return(data.frame(
      beta_adv = as.numeric(beta_adv), 
      se_adv = as.numeric(se_adv), 
      beta_unknown = as.numeric(beta_unknown)
    ))
    
  }, error = function(e) {
    # Silently fail to NA if a specific iteration crashes
    return(data.frame(beta_adv = NA, se_adv = NA, beta_unknown = NA))
  })
}


# ##############################################################################
# MODULE 4: Monte Carlo Engine
# ##############################################################################

# --- Block 4.1: Execution Loop ---
run_monte_carlo <- function(cohort, n_iter = 100) {
  cat("\n>>> MODULE 4, BLOCK 4.1: Running Monte Carlo Loop\n")
  cat(sprintf("Iterations per mechanism: %d\n", n_iter))
  
  # FIX: Set seed to ensure 100% reproducibility of the exact table
  set.seed(123)
  
  results <- list()
  
  for (mech in c("MCAR", "MAR", "MNAR")) {
    cat("\n  - Simulating mechanism:", mech, "\n")
    
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    
    for (i in 1:n_iter) {
      setTxtProgressBar(pb, i)
      
      if (mech == "MCAR") sim_data <- gen_mcar(cohort)
      if (mech == "MAR")  sim_data <- gen_mar(cohort)
      if (mech == "MNAR") sim_data <- gen_mnar(cohort)
      
      for (strat in c("S1", "S2", "S3", "S4")) {
        res <- evaluate_strategy(sim_data, strat)
        
        results[[length(results) + 1]] <- data.frame(
          Mechanism = mech,
          Iteration = i,
          Strategy = strat,
          Beta_Adv = res$beta_adv,
          SE_Adv = res$se_adv,
          Beta_Unknown = res$beta_unknown
        )
      }
    }
    close(pb)
  }
  
  results_df <- do.call(rbind, results)
  return(results_df)
}


# ##############################################################################
# MODULE 5: Summarization & Export
# ##############################################################################

# --- Block 5.1: Calculate Bias, MSE, Coverage ---
summarize_mc_results <- function(results_df, true_beta_adv, true_beta_unknown) {
  cat("\n>>> MODULE 5, BLOCK 5.1: Summarizing Results\n")
  
  summary_df <- results_df %>%
    group_by(Mechanism, Strategy) %>%
    summarize(
      Bias_Adv = mean(Beta_Adv - true_beta_adv, na.rm = TRUE),
      MSE_Adv = mean((Beta_Adv - true_beta_adv)^2, na.rm = TRUE),
      Coverage_Adv = mean((Beta_Adv - 1.96 * SE_Adv <= true_beta_adv) & 
                            (Beta_Adv + 1.96 * SE_Adv >= true_beta_adv), na.rm = TRUE),
      Bias_Unknown = mean(Beta_Unknown - true_beta_unknown, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(as.data.frame(summary_df))
  return(summary_df)
}

# --- Block 5.2: Save Results ---
save_results <- function(raw_results, summary_df) {
  cat("\n>>> MODULE 5, BLOCK 5.2: Saving Results\n")
  
  saveRDS(list(
    raw_results = raw_results,
    summary = summary_df
  ), "phase3a_sim_monte_carlo_results.rds")
  
  cat("Saved phase3a_sim_monte_carlo_results.rds successfully.\n")
}


# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_simulation_pipeline <- function() {
  cat("========================================================\n")
  cat("Starting Phase 3A Monte Carlo Simulation\n")
  cat("========================================================\n")
  
  init_environment()
  setup <- load_ground_truth()
  
  raw_results <- run_monte_carlo(setup$cohort, n_iter = 100)
  
  summary_df <- summarize_mc_results(raw_results, setup$true_beta_adv, setup$true_beta_unknown)
  
  save_results(raw_results, summary_df)
  
  cat("\n========================================================\n")
  cat("Simulation Pipeline Complete!\n")
  cat("========================================================\n")
}

# Execute the Master Pipeline unless a sourcing script has opted out.
# Scripts 10 and 11 source this file only for its missingness generators and
# evaluate_strategy(); they set options(phase3a.mc_execute = FALSE) before
# sourcing to prevent the full Monte Carlo from re-running.
if (isTRUE(getOption("phase3a.mc_execute", default = TRUE))) {
  run_simulation_pipeline()
}