# ==============================================================================
# Script 01: Phase 3A Empirical Arm - Expression Modeling (A vs B)
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
#
# Purpose:  Fit the high-dimensional prognostic models (LASSO, Elastic Net, 
#           Ridge, RSF, Stepwise Cox) on the Expression cohort using Strategy A 
#           (Complete-Case) and Strategy B (Missing-Retained). Extract selected 
#           features and calculate performance metrics.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Setup & Data Loading
#   MODULE 2: Strategy A Pipeline (Complete-Case)
#   MODULE 3: Strategy B Pipeline (Missing-Retained)
#   MODULE 4: Performance Evaluation & Feature Extraction
#   MODULE 5: Serialization & Export
#
# Save as:  01_empirical_expression.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/Phase_3A/
#
# Input:    phase3a_base.rds
# Output:   phase3a_empirical_expr.rds
# ==============================================================================

# ##############################################################################
# MODULE 1: Setup & Data Loading
# ##############################################################################

# --- Block 1.1: Load Data & Packages ---
init_phase3a_expr <- function() {
  cat("\n>>> MODULE 1, BLOCK 1.1: Loading Data and Packages\n")
  
  if (!file.exists("phase3a_base.rds")) stop("Run 00_base_setup_phase3a.R first!")
  base_data <- readRDS("phase3a_base.rds")
  
  # Libraries
  library(glmnet)
  library(survival)
  library(randomForestSRC)
  library(timeROC)
  
  return(base_data)
}

# --- Block 1.2: Core Modeling Functions ---
get_modeling_functions <- function() {
  cat(">>> MODULE 1, BLOCK 1.2: Loading Modeling Functions\n")
  
  # 1. Penalized Cox (LASSO, EN, Ridge)
  fit_penalized <- function(train_df, test_df, genes, alpha_val) {
    x_train <- as.matrix(train_df[, genes])
    y_train <- Surv(train_df$os_time, train_df$os_event)
    
    x_test  <- as.matrix(test_df[, genes])
    y_test  <- Surv(test_df$os_time, test_df$os_event)
    
    set.seed(123)
    cv_fit <- cv.glmnet(x_train, y_train, family = "cox", alpha = alpha_val, 
                        nfolds = 10, type.measure = "C")
    
    lp_test <- as.numeric(predict(cv_fit, newx = x_test, s = "lambda.min"))
    
    # Inverted C-index to correct package orientation
    cindex <- 1 - survival::concordance(y_test ~ lp_test)$concordance
    
    coefs <- coef(cv_fit, s = "lambda.min")
    active_genes <- rownames(coefs)[which(coefs != 0)]
    
    return(list(
      model = cv_fit,
      cindex = cindex,
      active_genes = active_genes,
      n_active = length(active_genes),
      lambdas = cv_fit$lambda.min
    ))
  }
  
  # 2. Random Survival Forest
  fit_rsf <- function(train_df, test_df, genes) {
    df_train <- train_df[, c("os_time", "os_event", genes)]
    df_test  <- test_df[, c("os_time", "os_event", genes)]
    
    set.seed(123)
    rsf_fit <- rfsrc(Surv(os_time, os_event) ~ ., data = df_train, 
                     ntree = 500, nodesize = 10, mtry = 5, importance = TRUE)
    
    preds <- predict(rsf_fit, newdata = df_test)
    cindex <- preds$err.rate[length(preds$err.rate)] 
    
    vimp <- sort(rsf_fit$importance, decreasing = TRUE)
    top_genes <- names(vimp)[1:5]
    
    return(list(model = rsf_fit, cindex = 1 - cindex, top_genes = top_genes))
  }
  
  # 3. Stepwise Cox (Clinical Base + Genes)
  fit_stepwise <- function(train_df, test_df, genes) {
    base_formula <- "Surv(os_time, os_event) ~ age + gender + histology + stage_group"
    
    # Wrap genes in backticks to handle hyphens (e.g., MIR124-1)
    safe_genes <- paste0("`", genes, "`")
    full_formula <- as.formula(paste(base_formula, "+", paste(safe_genes, collapse = " + ")))
    
    fit_base <- coxph(as.formula(base_formula), data = train_df)
    fit_full <- coxph(full_formula, data = train_df)
    
    lower_scope <- as.formula("~ age + gender + histology + stage_group")
    upper_scope <- as.formula(paste("~", paste(c("age", "gender", "histology", "stage_group", safe_genes), collapse = " + ")))
    
    fit_step <- step(fit_base, scope = list(lower = lower_scope, upper = upper_scope), 
                     direction = "both", trace = 0)
    
    lp_test <- as.numeric(predict(fit_step, newdata = test_df, type = "lp"))
    cindex <- 1 - survival::concordance(Surv(test_df$os_time, test_df$os_event) ~ lp_test)$concordance
    
    active_vars <- names(coef(fit_step))
    
    return(list(model = fit_step, cindex = cindex, active_vars = active_vars))
  }
  
  return(list(fit_penalized = fit_penalized, fit_rsf = fit_rsf, fit_stepwise = fit_stepwise))
}


# ##############################################################################
# MODULE 2: Strategy A Pipeline (Complete-Case, N=744)
# ##############################################################################

run_strategy_A_expr <- function(base_data, model_funcs) {
  cat("\n>>> MODULE 2: Running Strategy A (Complete-Case)\n")
  
  train <- base_data$train_A_expr
  test  <- base_data$test_A_expr
  genes <- base_data$present_genes
  
  results <- list()
  
  # --- Block 2.2: Penalized Models ---
  cat("  - Block 2.2: Penalized Models\n")
  results$lasso <- model_funcs$fit_penalized(train, test, genes, alpha_val = 1.0)
  results$enet  <- model_funcs$fit_penalized(train, test, genes, alpha_val = 0.5)
  results$ridge <- model_funcs$fit_penalized(train, test, genes, alpha_val = 0.0)
  
  # --- Block 2.3: RSF & Stepwise ---
  cat("  - Block 2.3: RSF and Stepwise\n")
  results$rsf <- model_funcs$fit_rsf(train, test, genes)
  results$stepwise <- model_funcs$fit_stepwise(train, test, genes)
  
  # Clinical Only Baseline
  clin_fit <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_group, data = train)
  lp_clin <- as.numeric(predict(clin_fit, newdata = test, type = "lp"))
  
  results$clinical <- list(
    model = clin_fit, 
    cindex = 1 - survival::concordance(Surv(test$os_time, test$os_event) ~ lp_clin)$concordance
  )
  
  return(results)
}


# ##############################################################################
# MODULE 3: Strategy B Pipeline (Missing-Retained, N=867)
# ##############################################################################

run_strategy_B_expr <- function(base_data, model_funcs) {
  cat("\n>>> MODULE 3: Running Strategy B (Missing-Retained)\n")
  
  train <- base_data$train_B_expr
  test  <- base_data$test_B_expr
  genes <- base_data$present_genes
  
  results <- list()
  
  # --- Block 3.2: Penalized Models ---
  cat("  - Block 3.2: Penalized Models\n")
  results$lasso <- model_funcs$fit_penalized(train, test, genes, alpha_val = 1.0)
  results$enet  <- model_funcs$fit_penalized(train, test, genes, alpha_val = 0.5)
  results$ridge <- model_funcs$fit_penalized(train, test, genes, alpha_val = 0.0)
  
  # --- Block 3.3: RSF & Stepwise ---
  cat("  - Block 3.3: RSF and Stepwise\n")
  results$rsf <- model_funcs$fit_rsf(train, test, genes)
  results$stepwise <- model_funcs$fit_stepwise(train, test, genes)
  
  # Clinical Only Baseline
  clin_fit <- coxph(Surv(os_time, os_event) ~ age + gender + histology + stage_group, data = train)
  lp_clin <- as.numeric(predict(clin_fit, newdata = test, type = "lp"))
  
  results$clinical <- list(
    model = clin_fit, 
    cindex = 1 - survival::concordance(Surv(test$os_time, test$os_event) ~ lp_clin)$concordance
  )
  
  return(results)
}


# ##############################################################################
# MODULE 4: Performance Evaluation & Feature Extraction
# ##############################################################################

summarize_results <- function(results_A, results_B) {
  cat("\n>>> MODULE 4: Summarizing Results\n")
  
  # Create a comparison table for C-indices
  cindex_table <- data.frame(
    Model = c("Clinical Only", "LASSO", "Elastic Net", "Ridge", "RSF", "Stepwise Cox"),
    Strategy_A_Cindex = c(
      results_A$clinical$cindex,
      results_A$lasso$cindex,
      results_A$enet$cindex,
      results_A$ridge$cindex,
      results_A$rsf$cindex,
      results_A$stepwise$cindex
    ),
    Strategy_B_Cindex = c(
      results_B$clinical$cindex,
      results_B$lasso$cindex,
      results_B$enet$cindex,
      results_B$ridge$cindex,
      results_B$rsf$cindex,
      results_B$stepwise$cindex
    )
  )
  
  # Feature Selection Comparison (LASSO)
  lasso_compare <- list(
    Strategy_A_Genes = results_A$lasso$active_genes,
    Strategy_B_Genes = results_B$lasso$active_genes,
    Intersection = intersect(results_A$lasso$active_genes, results_B$lasso$active_genes),
    A_Only = setdiff(results_A$lasso$active_genes, results_B$lasso$active_genes),
    B_Only = setdiff(results_B$lasso$active_genes, results_A$lasso$active_genes)
  )
  
  cat("\n--- C-Index Comparison ---\n")
  print(cindex_table, row.names = FALSE)
  
  cat("\n--- LASSO Feature Selection Comparison ---\n")
  cat("Strategy A (Complete-Case) selected:", length(lasso_compare$Strategy_A_Genes), "genes\n")
  cat("Strategy B (Missing-Retained) selected:", length(lasso_compare$Strategy_B_Genes), "genes\n")
  cat("Intersection:", paste(lasso_compare$Intersection, collapse = ", "), "\n")
  cat("Dropped by A (B only):", paste(lasso_compare$B_Only, collapse = ", "), "\n")
  cat("Spurious in A (A only):", paste(lasso_compare$A_Only, collapse = ", "), "\n")
  
  return(list(cindex_table = cindex_table, lasso_compare = lasso_compare))
}


# ##############################################################################
# MODULE 5: Serialization & Export
# ##############################################################################

save_empirical_results <- function(results_A, results_B, summary_res) {
  cat("\n>>> MODULE 5: Saving Results\n")
  
  saveRDS(list(
    results_A = results_A,
    results_B = results_B,
    cindex_table = summary_res$cindex_table,
    lasso_compare = summary_res$lasso_compare
  ), "phase3a_empirical_expr.rds")
  
  cat("Saved phase3a_empirical_expr.rds successfully.\n")
}


# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_empirical_expression_pipeline <- function() {
  cat("========================================================\n")
  cat("Starting Phase 3A Empirical Expression Pipeline\n")
  cat("========================================================\n")
  
  # Module 1
  base_data <- init_phase3a_expr()
  model_funcs <- get_modeling_functions()
  
  # Module 2
  res_A <- run_strategy_A_expr(base_data, model_funcs)
  
  # Module 3
  res_B <- run_strategy_B_expr(base_data, model_funcs)
  
  # Module 4
  summary_res <- summarize_results(res_A, res_B)
  
  # Module 5
  save_empirical_results(res_A, res_B, summary_res)
  
  cat("\n========================================================\n")
  cat("Phase 3A Empirical Expression Pipeline Complete!\n")
  cat("========================================================\n")
}

# Execute
run_empirical_expression_pipeline()