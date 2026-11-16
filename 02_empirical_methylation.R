# ==============================================================================
# Script 02: Phase 3A Empirical Arm - Methylation Modeling (A vs B)
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
#
# Purpose:  Fit the high-dimensional prognostic models on the Methylation 
#           cohort (N=701) using Strategy A (Complete-Case) and Strategy B 
#           (Missing-Retained). Evaluates Gene-Level models.
#           Calculates validation C-index and HR per IQR (Full Cohort).
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Setup & Data Loading
#   MODULE 2: Strategy A Pipeline (Complete-Case)
#   MODULE 3: Strategy B Pipeline (Missing-Retained)
#   MODULE 4: Performance Evaluation & Feature Extraction
#   MODULE 5: Serialization & Export
#
# Input:    phase3a_base.rds
# Output:   phase3a_empirical_meth.rds
# ==============================================================================


# ##############################################################################
# MODULE 1: Setup & Data Loading
# ##############################################################################

# --- Block 1.1: Load Data & Packages ---
init_meth_pipeline <- function() {
  cat("\n>>> MODULE 1, BLOCK 1.1: Loading Data and Packages\n")
  
  if (!file.exists("phase3a_base.rds")) stop("Run 00_base_setup_phase3a.R first!")
  base_data <- readRDS("phase3a_base.rds")
  
  library(glmnet)
  library(survival)
  
  return(base_data)
}

# --- Block 1.2: Methylation Modeling Functions ---
get_meth_functions <- function() {
  cat(">>> MODULE 1, BLOCK 1.2: Loading Modeling Functions\n")
  
  # Function to fit LASSO to methylation data and return C-index + HR per IQR
  fit_meth_lasso <- function(train_df, test_df, features) {
    features <- intersect(features, colnames(train_df))
    features <- intersect(features, colnames(test_df))
    
    x_train <- as.matrix(train_df[, features])
    y_train <- Surv(train_df$os_time, train_df$os_event)
    
    x_test  <- as.matrix(test_df[, features])
    y_test  <- Surv(test_df$os_time, test_df$os_event)
    
    # Fit LASSO
    set.seed(123)
    cv_fit <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1.0, 
                        nfolds = 10, type.measure = "C")
    
    # Extract coefficients safely
    coefs <- coef(cv_fit, s = "lambda.min")
    active_idx <- which(coefs != 0)
    active_features <- rownames(coefs)[active_idx]
    
    if (length(active_features) == 0) {
      return(list(cindex = 0.5, hr_iqr = NA, n_active = 0, active_features = character(0)))
    }
    
    # Calculate Risk Score (FIX: using as.matrix to avoid sparse matrix issues)
    beta <- as.numeric(as.matrix(coefs)[active_idx, 1])
    rs_train <- as.numeric(x_train[, active_features, drop = FALSE] %*% beta)
    iqr_val <- IQR(rs_train, na.rm = TRUE) # FIX: Added na.rm
    rs_test <- as.numeric(x_test[, active_features, drop = FALSE] %*% beta)
    
    # Validation C-index (Inverted)
    cindex <- 1 - survival::concordance(y_test ~ rs_test)$concordance
    
    # Multivariable Cox on FULL COHORT for HR per IQR
    full_eval <- rbind(train_df, test_df)
    x_full <- as.matrix(full_eval[, active_features, drop = FALSE])
    rs_full <- as.numeric(x_full %*% beta)
    
    full_eval$risk_score <- rs_full
    
    cox_adj <- coxph(Surv(os_time, os_event) ~ risk_score + age + gender + histology + stage_group, data = full_eval)
    
    beta_adj <- coef(cox_adj)["risk_score"]
    se_adj   <- sqrt(vcov(cox_adj)["risk_score", "risk_score"])
    
    hr_iqr <- exp(beta_adj * iqr_val)
    hr_lower <- exp((beta_adj - 1.96 * se_adj) * iqr_val)
    hr_upper <- exp((beta_adj + 1.96 * se_adj) * iqr_val)
    p_val <- summary(cox_adj)$coefficients["risk_score", "Pr(>|z|)"]
    
    return(list(
      model = cv_fit,
      cindex = cindex,
      hr_iqr = hr_iqr,
      hr_ci = c(hr_lower, hr_upper),
      p_value = p_val,
      n_active = length(active_features),
      active_features = active_features,
      iqr_value = iqr_val
    ))
  }
  
  return(list(fit_meth_lasso = fit_meth_lasso))
}


# ##############################################################################
# MODULE 2: Strategy A Pipeline (Complete-Case)
# ##############################################################################

run_strategy_A_meth <- function(base_data, meth_funcs) {
  cat("\n>>> MODULE 2: Running Strategy A Methylation (Complete-Case)\n")
  
  train <- base_data$train_A_meth
  test  <- base_data$test_A_meth
  
  # Gene-level features (24 genes)
  gene_features <- base_data$anno_prom$UCSC_RefGene_Name
  gene_features <- unique(unlist(strsplit(gene_features, ";")))
  gene_features <- intersect(gene_features, colnames(train))
  
  cat("  - Fitting Gene-Level LASSO\n")
  res_gene <- meth_funcs$fit_meth_lasso(train, test, gene_features)
  
  # Clinical Baseline (Inverted C-index)
  clin_fit <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_group, data = train)
  lp_clin <- as.numeric(predict(clin_fit, newdata = test, type = "lp"))
  cindex_clin <- 1 - survival::concordance(Surv(test$os_time, test$os_event) ~ lp_clin)$concordance
  
  return(list(
    gene_level = res_gene,
    clinical_cindex = cindex_clin
  ))
}


# ##############################################################################
# MODULE 3: Strategy B Pipeline (Missing-Retained)
# ##############################################################################

run_strategy_B_meth <- function(base_data, meth_funcs) {
  cat("\n>>> MODULE 3: Running Strategy B Methylation (Missing-Retained)\n")
  
  train <- base_data$train_B_meth
  test  <- base_data$test_B_meth
  
  # Gene-level features (24 genes)
  gene_features <- base_data$anno_prom$UCSC_RefGene_Name
  gene_features <- unique(unlist(strsplit(gene_features, ";")))
  gene_features <- intersect(gene_features, colnames(train))
  
  cat("  - Fitting Gene-Level LASSO\n")
  res_gene <- meth_funcs$fit_meth_lasso(train, test, gene_features)
  
  # Clinical Baseline (Inverted C-index)
  clin_fit <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_group, data = train)
  lp_clin <- as.numeric(predict(clin_fit, newdata = test, type = "lp"))
  cindex_clin <- 1 - survival::concordance(Surv(test$os_time, test$os_event) ~ lp_clin)$concordance
  
  return(list(
    gene_level = res_gene,
    clinical_cindex = cindex_clin
  ))
}


# ##############################################################################
# MODULE 4: Performance Evaluation
# ##############################################################################

summarize_meth_results <- function(res_A, res_B) {
  cat("\n>>> MODULE 4: Summarizing Methylation Results\n")
  
  cat("\n--- C-Index Comparison (Gene-Level Methylation) ---\n")
  cat(sprintf("Strategy A (Complete-Case): Validation C-index = %.3f (Clinical Baseline = %.3f)\n", 
              res_A$gene_level$cindex, res_A$clinical_cindex))
  cat(sprintf("Strategy B (Missing-Retained): Validation C-index = %.3f (Clinical Baseline = %.3f)\n", 
              res_B$gene_level$cindex, res_B$clinical_cindex))
  
  cat("\n--- Hazard Ratio per IQR (Adjusted for Clinical, Full Cohort) ---\n")
  
  if(!is.na(res_A$gene_level$hr_iqr)) {
    cat(sprintf("Strategy A: HR = %.2f (95%% CI: %.2f-%.2f), p = %.3f\n", 
                res_A$gene_level$hr_iqr, res_A$gene_level$hr_ci[1], res_A$gene_level$hr_ci[2], res_A$gene_level$p_value))
  } else {
    cat("Strategy A: No features selected.\n")
  }
  
  if(!is.na(res_B$gene_level$hr_iqr)) {
    cat(sprintf("Strategy B: HR = %.2f (95%% CI: %.2f-%.2f), p = %.3f\n", 
                res_B$gene_level$hr_iqr, res_B$gene_level$hr_ci[1], res_B$gene_level$hr_ci[2], res_B$gene_level$p_value))
  } else {
    cat("Strategy B: No features selected.\n")
  }
  
  cat("\n--- LASSO Feature Selection Comparison (Gene-Level) ---\n")
  a_genes <- res_A$gene_level$active_features
  b_genes <- res_B$gene_level$active_features
  
  cat("Strategy A selected:", length(a_genes), "genes\n")
  cat("Strategy B selected:", length(b_genes), "genes\n")
  cat("Intersection:", paste(intersect(a_genes, b_genes), collapse = ", "), "\n")
  cat("Unique to Strategy A:", paste(setdiff(a_genes, b_genes), collapse = ", "), "\n")
  cat("Unique to Strategy B:", paste(setdiff(b_genes, a_genes), collapse = ", "), "\n")
}


# ##############################################################################
# MODULE 5: Serialization & Export
# ##############################################################################

save_meth_results <- function(res_A, res_B) {
  cat("\n>>> MODULE 5: Saving Methylation Results\n")
  
  saveRDS(list(
    results_A = res_A,
    results_B = res_B
  ), "phase3a_empirical_meth.rds")
  
  cat("Saved phase3a_empirical_meth.rds successfully.\n")
}


# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_empirical_methylation_pipeline <- function() {
  cat("========================================================\n")
  cat("Starting Phase 3A Empirical Methylation Pipeline\n")
  cat("========================================================\n")
  
  base_data <- init_meth_pipeline()
  meth_funcs <- get_meth_functions()
  
  res_A <- run_strategy_A_meth(base_data, meth_funcs)
  res_B <- run_strategy_B_meth(base_data, meth_funcs)
  
  summarize_meth_results(res_A, res_B)
  save_meth_results(res_A, res_B)
  
  cat("\n========================================================\n")
  cat("Phase 3A Empirical Methylation Pipeline Complete!\n")
  cat("========================================================\n")
}

# Execute
run_empirical_methylation_pipeline()