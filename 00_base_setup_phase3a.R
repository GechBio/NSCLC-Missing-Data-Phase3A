# ==============================================================================
# Script 00: Phase 3A Shared Setup, Utilities, and Dual-Cohort Assembly
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Build the analytic NSCLC cohort from TCGA RNA-seq and clinical
#           data; construct promoter-region gene-averaged methylation matrices;
#           define shared utilities; define Strategy A (Complete-Case) and 
#           Strategy B (Missing-Retained) cohorts; perform stratified 70/30 
#           splits for both strategies; save all base objects to disk.
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#       - Block 1.1: Path Definition & Directory Creation
#       - Block 1.2: Package Installation & Loading
#   MODULE 2: Shared Statistical & Parsing Utilities
#       - Block 2.1: Modeling Utilities
#       - Block 2.2: Parsing Utilities
#   MODULE 3: Expression Cohort Construction
#       - Block 3.1 & 3.2: Load Expression Data & Format Clinical Covariates
#       - Block 3.3: Merge & Filter (Master Expression N=867)
#   MODULE 4: Methylation Cohort Construction
#       - Block 4.1 & 4.2: Load, Map, and Filter Methylation Probes
#       - Block 4.3: Merge & Filter (Master Methylation N=701)
#   MODULE 5: Strategy Definition & Data Splitting
#       - Block 5.1 & 5.2: Define A/B and Split
#   MODULE 6: Serialization & Export
#       - Block 6.1 & 6.2: Package and Save
#
# Save as:  00_base_setup_phase3a.R
# Save in:  ~/TCGA_Local_Storage/TCGA NSCLC Project/Phase_3A/
#
# Output:   phase3a_base.rds
# Seeds:    split = 123
# ==============================================================================


# ##############################################################################
# MODULE 1: Environment & Dependency Management
# ##############################################################################

# --- Block 1.1: Path Definition & Directory Creation (portable) ---
init_paths <- function() {
  cat("\n>>> MODULE 1, BLOCK 1.1: Setting Paths\n")
  
  resolve_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- grep("^--file=", args, value = TRUE)
    if (length(f) > 0) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
    if (!is.null(sys.frame(1)$ofile)) {
      return(dirname(normalizePath(sys.frame(1)$ofile)))
    }
    stop("Cannot determine script location. Run via Rscript or set PROJECT_ROOT.")
  }
  
  script_dir   <- resolve_script_dir()
  project_root <- Sys.getenv("PROJECT_ROOT",
                             unset = normalizePath(file.path(script_dir, "..")))
  
  paths <- list(
    parent_dir  = project_root,
    phase3a_dir = file.path(project_root, "Phase_3A"),
    meth_dir    = Sys.getenv("METH_DIR",
                             unset = file.path(project_root, "nsclc_raw_data"))
  )
  paths$sample_sheet <- file.path(paths$parent_dir, "NSCLC_GDC_Sample_Sheet.tsv")
  paths$expr_rds     <- file.path(paths$parent_dir, "master_nsclc_matrix.rds")
  paths$surv_rds     <- file.path(paths$parent_dir, "clinical_final.rds")
  paths$cov_rds      <- file.path(paths$parent_dir, "clinical_covariates.rds")
  
  if (!dir.exists(paths$phase3a_dir)) {
    dir.create(paths$phase3a_dir, recursive = TRUE)
    cat("Created Phase 3A directory:", paths$phase3a_dir, "\n")
  }
  setwd(paths$phase3a_dir)
  
  # Validate the four parent-pipeline files
  required <- c(paths$expr_rds, paths$surv_rds, paths$cov_rds, paths$sample_sheet)
  missing  <- required[!file.exists(required)]
  if (length(missing) > 0) {
    stop("Missing parent-pipeline files:\n  ",
         paste(missing, collapse = "\n  "),
         "\nRun the Paper 2 pipeline (00_base_setup.R) first, or set PROJECT_ROOT.")
  }
  
  # Validate the methylation directory
  meth_files <- list.files(paths$meth_dir,
                           pattern = "\\.methylation_array\\.sesame\\.level3betas\\.txt$",
                           recursive = TRUE, full.names = TRUE)
  if (length(meth_files) == 0) {
    stop("No methylation files found in METH_DIR (", paths$meth_dir, ").\n",
         "Set METH_DIR to the directory containing ",
         "*.methylation_array.sesame.level3betas.txt files, or copy them there ",
         "before running this script.")
  }
  cat("Methylation files found:", length(meth_files), "\n")
  
  paths
}

# --- Block 1.2: Package Installation & Loading ---
load_packages <- function() {
  cat(">>> MODULE 1, BLOCK 1.2: Loading Packages\n")
  
  pkgs <- c("glmnet","survival","survminer","timeROC","rms","dcurves","caret",
            "dplyr","tidyr","ggplot2","patchwork","pheatmap","data.table","tibble",
            "randomForestSRC","purrr","magick",
            "IlluminaHumanMethylation450kanno.ilmn12.hg19")
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (p == "IlluminaHumanMethylation450kanno.ilmn12.hg19") {
        BiocManager::install(p, update = FALSE, ask = FALSE)
      } else install.packages(p)
    }
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}


# ##############################################################################
# MODULE 2: Shared Statistical & Parsing Utilities
# ##############################################################################

# --- Block 2.1: Modeling Utilities ---
get_model_utils <- function() {
  cat(">>> MODULE 2, BLOCK 2.1: Loading Modeling Utilities\n")
  
  calc_risk <- function(data, genes, coefs) {
    genes <- intersect(genes, colnames(data))
    if (length(genes) == 0) return(rep(NA_real_, nrow(data)))
    cf <- as.numeric(coefs[genes])
    as.numeric(as.matrix(data[, genes, drop = FALSE]) %*% cf)
  }
  
  safe_cindex <- function(train_fit, newdata, time, event) {
    lp <- tryCatch(predict(train_fit, newdata = newdata, type = "lp"),
                   error = function(e) rep(NA_real_, nrow(newdata)))
    if (all(is.na(lp))) return(NA_real_)
    cc <- tryCatch(
      unname(summary(survival::coxph(survival::Surv(time, event) ~ lp))$concordance[1]),
      error = function(e) NA_real_)
    as.numeric(cc)
  }
  
  rescale_hr <- function(beta, se, scale_factor, scale_label = "per unit") {
    if (!is.numeric(beta) || length(beta) != 1 || is.na(beta)) stop("rescale_hr: beta must be a single non-NA numeric value")
    if (!is.numeric(se) || length(se) != 1 || is.na(se)) stop("rescale_hr: se must be a single non-NA numeric value")
    if (!is.numeric(scale_factor) || length(scale_factor) != 1 || is.na(scale_factor) || scale_factor == 0) stop("rescale_hr: scale_factor must be a single non-zero numeric value")
    
    hr    <- exp(beta * scale_factor)
    ci_lo <- exp((beta - 1.96 * se) * scale_factor)
    ci_hi <- exp((beta + 1.96 * se) * scale_factor)
    p     <- 2 * pnorm(-abs(beta / se))
    
    data.frame(Scale = scale_label, Factor = scale_factor, HR = hr, 
               CI_lower = ci_lo, CI_upper = ci_hi, P_value = p, stringsAsFactors = FALSE)
  }
  
  return(list(calc_risk = calc_risk, safe_cindex = safe_cindex, rescale_hr = rescale_hr))
}

# --- Block 2.2: Parsing Utilities ---
get_parse_utils <- function() {
  cat(">>> MODULE 2, BLOCK 2.2: Loading Parsing Utilities\n")
  
  parse_stage <- function(x) {
    xu <- toupper(trimws(as.character(x)))
    out <- dplyr::case_when(
      grepl("IV",  xu) ~ "Advanced",
      grepl("III", xu) ~ "Advanced",
      grepl("II",  xu) ~ "Early",
      grepl("I",   xu) ~ "Early",
      TRUE ~ "Unknown"
    )
    factor(out, levels = c("Early","Advanced","Unknown"))
  }
  
  return(list(parse_stage = parse_stage))
}


# ##############################################################################
# MODULE 3: Expression Cohort Construction
# ##############################################################################

# --- Block 3.1 & 3.2: Load Expression Data & Format Clinical Covariates ---
load_and_format_expression <- function(paths, parse_stage) {
  cat(">>> MODULE 3, BLOCK 3.1 & 3.2: Loading and Formatting Expression/Clinical Data\n")
  
  master_nsclc_matrix <- readRDS(paths$expr_rds)
  clinical_final      <- readRDS(paths$surv_rds)
  clinical_covariates <- readRDS(paths$cov_rds)
  
  candidate_genes <- c(
    "APC","BRMS1","BTG2","CDKN2A","CHRNB4","DAPK1","EGLN2",
    "HOXA9","HDGFL3","KDM1A","KDM2A","LRRC3B","MGMT",
    "PCDH10","PITX2","RASSF1","RUNX3","SHOX2","SPARC",
    "TGFBI","WIF1","MIR34B","MIR124-1","MIR124-2","MIR124-3"
  )
  
  surv_endpoints <- data.frame(
    patient_id             = clinical_final$cases.submitter_id,
    vital_status           = clinical_final$demographic.vital_status,
    days_to_death          = suppressWarnings(as.numeric(clinical_final$demographic.days_to_death)),
    days_to_last_follow_up = suppressWarnings(as.numeric(clinical_final$diagnoses.days_to_last_follow_up)),
    stringsAsFactors = FALSE)
  surv_endpoints$os_time  <- ifelse(surv_endpoints$vital_status == "Dead",
                                    surv_endpoints$days_to_death,
                                    surv_endpoints$days_to_last_follow_up)
  surv_endpoints$os_event <- as.integer(surv_endpoints$vital_status == "Dead")
  surv_endpoints <- surv_endpoints[!is.na(surv_endpoints$os_time) & surv_endpoints$os_time >= 30, ]
  
  covariates_clean <- data.frame(
    patient_id = clinical_covariates$cases.submitter_id,
    age        = suppressWarnings(as.numeric(clinical_covariates$demographic.age_at_index)),
    gender     = tolower(clinical_covariates$demographic.gender),
    stage      = clinical_covariates$diagnoses.ajcc_pathologic_stage,
    project    = clinical_covariates$project.project_id,
    stringsAsFactors = FALSE)
  covariates_clean <- covariates_clean[grepl("TCGA-LUAD|TCGA-LUSC", covariates_clean$project), ]
  covariates_clean <- covariates_clean[!duplicated(covariates_clean$patient_id), ]
  covariates_clean$stage_group <- parse_stage(covariates_clean$stage)
  covariates_clean$histology   <- ifelse(grepl("LUAD", covariates_clean$project), "LUAD",
                                         ifelse(grepl("LUSC", covariates_clean$project), "LUSC", NA))
  covariates_clean$histology   <- factor(covariates_clean$histology, levels = c("LUAD","LUSC"))
  covariates_clean$gender      <- factor(covariates_clean$gender, levels = c("female","male"))
  
  return(list(
    matrix = master_nsclc_matrix, 
    surv = surv_endpoints, 
    covs = covariates_clean, 
    candidate_genes = candidate_genes
  ))
}

# --- Block 3.3: Merge & Filter (Master Expression N=867) ---
merge_expression_cohort <- function(expr_raw) {
  cat(">>> MODULE 3, BLOCK 3.3: Merging and Filtering Expression Cohort\n")
  
  present_genes <- intersect(expr_raw$candidate_genes, rownames(expr_raw$matrix))
  sample_ids    <- colnames(expr_raw$matrix)
  expr_t        <- as.data.frame(t(expr_raw$matrix[present_genes, , drop = FALSE]))
  expr_t$patient_id <- sample_ids
  
  cohort <- merge(expr_t, expr_raw$surv, by = "patient_id")
  cohort <- merge(cohort, expr_raw$covs, by = "patient_id", all.x = TRUE)
  cohort <- cohort[!is.na(cohort$stage_group) & !is.na(cohort$histology) & 
                     !is.na(cohort$age) & !is.na(cohort$gender), ]
  cohort <- cohort[order(cohort$patient_id, -cohort$os_time), ]
  cohort <- cohort[!duplicated(cohort$patient_id), ]
  
  return(list(
    cohort = cohort, 
    present_genes = present_genes,
    surv = expr_raw$surv,
    covs = expr_raw$covs,
    candidate_genes = expr_raw$candidate_genes
  ))
}


# ##############################################################################
# MODULE 4: Methylation Cohort Construction
# ##############################################################################

# --- Block 4.1 & 4.2: Load, Map, and Filter Methylation Probes ---
build_methylation_matrix <- function(paths, candidate_genes) {
  cat(">>> MODULE 4, BLOCK 4.1 & 4.2: Loading and Processing Methylation Matrices\n")
  
  meth_files <- list.files(paths$meth_dir, pattern = "\\.methylation_array\\.sesame\\.level3betas\\.txt$", recursive = TRUE, full.names = TRUE)
  if (length(meth_files) == 0) stop("No methylation files found.")
  
  sheet <- read.table(paths$sample_sheet, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sheet_meth <- sheet[sheet$Data.Type == "Methylation Beta Value", ]
  folder_uuids <- basename(dirname(meth_files))
  uuid_to_patient <- setNames(sheet_meth$Case.ID, sheet_meth$File.ID)
  
  first_meth <- data.table::fread(meth_files[1], data.table = FALSE)
  probe_ids  <- first_meth[[1]]
  meth_beta <- matrix(NA_real_, nrow = length(probe_ids), ncol = length(meth_files), dimnames = list(probe_ids, folder_uuids))
  for (i in seq_along(meth_files)) {
    dat <- data.table::fread(meth_files[i], data.table = FALSE)
    meth_beta[, i] <- dat[[2]]
  }
  
  patient_ids <- uuid_to_patient[folder_uuids]
  colnames(meth_beta) <- patient_ids
  meth_beta <- meth_beta[, !is.na(colnames(meth_beta)), drop = FALSE]
  unique_patients <- unique(colnames(meth_beta))
  meth_patient <- matrix(NA_real_, nrow = nrow(meth_beta), ncol = length(unique_patients), dimnames = list(rownames(meth_beta), unique_patients))
  
  for (p in unique_patients) {
    cols <- which(colnames(meth_beta) == p)
    meth_patient[, p] <- if (length(cols) > 1) rowMeans(meth_beta[, cols, drop = FALSE], na.rm = TRUE) else meth_beta[, cols]
  }
  
  anno    <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  anno_df <- as.data.frame(anno)
  promoter_groups <- c("TSS1500","TSS200","5'UTR","1stExon")
  keep_region <- vapply(anno_df$UCSC_RefGene_Group, function(g) { if (is.na(g)) return(FALSE); any(strsplit(g, ";")[[1]] %in% promoter_groups) }, logical(1))
  anno_prom <- anno_df[keep_region, c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group")]
  anno_prom <- anno_prom[anno_prom$Name %in% rownames(meth_patient), ]
  
  probe_to_genes <- lapply(strsplit(anno_prom$UCSC_RefGene_Name, ";"), function(x) { x <- x[!is.na(x)]; intersect(x, candidate_genes) })
  keep_match <- lengths(probe_to_genes) > 0
  anno_prom  <- anno_prom[keep_match, ]
  probe_to_genes <- probe_to_genes[keep_match]
  
  gene_meth <- matrix(NA_real_, nrow = length(candidate_genes), ncol = ncol(meth_patient), dimnames = list(candidate_genes, colnames(meth_patient)))
  for (g in candidate_genes) {
    idx <- which(vapply(probe_to_genes, function(x) g %in% x, logical(1)))
    probes <- intersect(anno_prom$Name[idx], rownames(meth_patient))
    if (length(probes) >= 1) {
      gene_meth[g, ] <- colMeans(meth_patient[probes, , drop = FALSE], na.rm = TRUE)
    }
  }
  gene_meth <- gene_meth[rowSums(!is.na(gene_meth)) > 0, , drop = FALSE]
  
  return(list(gene_meth = gene_meth, anno_prom = anno_prom, probe_to_genes = probe_to_genes))
}

# --- Block 4.3: Merge & Filter (Master Methylation N=701) ---
merge_methylation_cohort <- function(meth_matrix_data, expr_cohort_data) {
  cat(">>> MODULE 4, BLOCK 4.3: Merging Methylation Cohort\n")
  
  gene_meth <- meth_matrix_data$gene_meth
  cov_meta <- expr_cohort_data$covs[, c("patient_id","age","gender","stage_group","histology")]
  
  meth_t   <- as.data.frame(t(gene_meth))
  meth_t$patient_id <- colnames(gene_meth)
  
  meth_cohort <- merge(meth_t, expr_cohort_data$surv, by = "patient_id")
  meth_cohort <- merge(meth_cohort, cov_meta, by = "patient_id", all.x = TRUE)
  meth_cohort <- meth_cohort[complete.cases(meth_cohort[, c(rownames(gene_meth), "os_time", "os_event", "age","gender","stage_group","histology")]), ]
  
  return(meth_cohort)
}


# ##############################################################################
# MODULE 5: Strategy Definition & Data Splitting
# ##############################################################################

# --- Block 5.1 & 5.2: Define A/B and Split ---
create_strategy_splits <- function(cohort_expr, cohort_meth) {
  cat(">>> MODULE 5, BLOCK 5.1 & 5.2: Defining Strategies and Splitting Data\n")
  
  # Define Strategies
  cohort_B_expr <- cohort_expr
  cohort_A_expr <- cohort_expr[cohort_expr$stage_group != "Unknown", ]
  cohort_A_expr$stage_group <- droplevels(cohort_A_expr$stage_group) 
  
  meth_cohort_B <- cohort_meth
  # BUG FIX APPLIED HERE: changed 'meth_cohort' to 'cohort_meth'
  meth_cohort_A <- cohort_meth[cohort_meth$stage_group != "Unknown", ]
  meth_cohort_A$stage_group <- droplevels(meth_cohort_A$stage_group)
  
  # Stratified Splits (Seed = 123)
  set.seed(123)
  
  strata_A_expr <- interaction(cohort_A_expr$os_event, cohort_A_expr$stage_group, cohort_A_expr$histology, drop = TRUE)
  idx_A_expr <- caret::createDataPartition(strata_A_expr, p = 0.7, list = FALSE)
  
  strata_B_expr <- interaction(cohort_B_expr$os_event, cohort_B_expr$stage_group, cohort_B_expr$histology, drop = TRUE)
  idx_B_expr <- caret::createDataPartition(strata_B_expr, p = 0.7, list = FALSE)
  
  strata_A_meth <- interaction(meth_cohort_A$os_event, meth_cohort_A$stage_group, meth_cohort_A$histology, drop = TRUE)
  idx_A_meth <- caret::createDataPartition(strata_A_meth, p = 0.7, list = FALSE)
  
  strata_B_meth <- interaction(meth_cohort_B$os_event, meth_cohort_B$stage_group, meth_cohort_B$histology, drop = TRUE)
  idx_B_meth <- caret::createDataPartition(strata_B_meth, p = 0.7, list = FALSE)
  
  return(list(
    cohort_A_expr = cohort_A_expr, cohort_B_expr = cohort_B_expr,
    train_A_expr = cohort_A_expr[idx_A_expr, ], test_A_expr = cohort_A_expr[-idx_A_expr, ],
    train_B_expr = cohort_B_expr[idx_B_expr, ], test_B_expr = cohort_B_expr[-idx_B_expr, ],
    meth_cohort_A = meth_cohort_A, meth_cohort_B = meth_cohort_B,
    train_A_meth = meth_cohort_A[idx_A_meth, ], test_A_meth = meth_cohort_A[-idx_A_meth, ],
    train_B_meth = meth_cohort_B[idx_B_meth, ], test_B_meth = meth_cohort_B[-idx_B_meth, ]
  ))
}


# ##############################################################################
# MODULE 6: Serialization & Export
# ##############################################################################

# --- Block 6.1 & 6.2: Package and Save ---
save_phase3a_base <- function(splits, expr_merged, meth_matrix, meth_cohort) {
  cat(">>> MODULE 6, BLOCK 6.1 & 6.2: Serializing Outputs\n")
  
  output_list <- list(
    # Cohorts
    cohort_A_expr = splits$cohort_A_expr,
    cohort_B_expr = splits$cohort_B_expr,
    meth_cohort_A = splits$meth_cohort_A,
    meth_cohort_B = splits$meth_cohort_B,
    
    # Splits
    train_A_expr = splits$train_A_expr, test_A_expr = splits$test_A_expr,
    train_B_expr = splits$train_B_expr, test_B_expr = splits$test_B_expr,
    train_A_meth = splits$train_A_meth, test_A_meth = splits$test_A_meth,
    train_B_meth = splits$train_B_meth, test_B_meth = splits$test_B_meth,
    
    # Metadata
    gene_meth = meth_matrix$gene_meth,
    anno_prom = meth_matrix$anno_prom,
    probe_to_genes = meth_matrix$probe_to_genes,
    candidate_genes = expr_merged$candidate_genes,
    present_genes = expr_merged$present_genes,
    covariates_clean = expr_merged$covs,
    surv_endpoints = expr_merged$surv
  )
  
  saveRDS(output_list, "phase3a_base.rds")
  cat("\nSaved phase3a_base.rds successfully in Phase_3A folder.\n")
  
  # Print Quick Summary
  cat("\n--- Final Cohort Sizes ---\n")
  cat(sprintf("Expression Strategy A (Complete-Case): Train=%d, Test=%d (Total=%d)\n", 
              nrow(splits$train_A_expr), nrow(splits$test_A_expr), nrow(splits$cohort_A_expr)))
  cat(sprintf("Expression Strategy B (Missing-Retained): Train=%d, Test=%d (Total=%d)\n", 
              nrow(splits$train_B_expr), nrow(splits$test_B_expr), nrow(splits$cohort_B_expr)))
  cat(sprintf("Methylation Strategy A (Complete-Case): Train=%d, Test=%d (Total=%d)\n", 
              nrow(splits$train_A_meth), nrow(splits$test_A_meth), nrow(splits$meth_cohort_A)))
  cat(sprintf("Methylation Strategy B (Missing-Retained): Train=%d, Test=%d (Total=%d)\n", 
              nrow(splits$train_B_meth), nrow(splits$test_B_meth), nrow(splits$meth_cohort_B)))
}


# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_phase3a_pipeline <- function() {
  cat("========================================================\n")
  cat("Starting Phase 3A Base Pipeline\n")
  cat("========================================================\n")
  
  # --- MODULE 1 ---
  paths <- init_paths()
  load_packages()
  
  # --- MODULE 2 ---
  model_utils <- get_model_utils()
  parse_utils <- get_parse_utils()
  
  # --- MODULE 3 ---
  expr_raw     <- load_and_format_expression(paths, parse_utils$parse_stage)
  expr_merged  <- merge_expression_cohort(expr_raw)
  
  # --- MODULE 4 ---
  meth_matrix  <- build_methylation_matrix(paths, expr_merged$candidate_genes)
  meth_cohort  <- merge_methylation_cohort(meth_matrix, expr_merged)
  
  # --- MODULE 5 ---
  splits       <- create_strategy_splits(expr_merged$cohort, meth_cohort)
  
  # --- MODULE 6 ---
  save_phase3a_base(splits, expr_merged, meth_matrix, meth_cohort)
  
  cat("\n========================================================\n")
  cat("Phase 3A Base Pipeline Complete!\n")
  cat("========================================================\n")
}

# Execute the master pipeline
run_phase3a_pipeline()