# ==============================================================================
# Script 03b: Generate Gene Selection Dot Plot 
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
# Contact:  derbew1953@gmail.com / d610113006@tmu.edu.tw
#
# Purpose:  Generate a presence/absence dot plot to visualize the intersections 
#           and unique sets of selected genes. 
#
# ARCHITECTURE: MODULE -> BLOCK DESIGN
#   MODULE 1: Environment & Dependency Management
#   MODULE 2: Data Preparation
#   MODULE 3: Plot Generation & Export
#
# Input:    phase3a_empirical_expr.rds
#           phase3a_empirical_meth.rds
# Output:   Figure3_Feature_Selection_DotPlot.pdf
# ==============================================================================

# ##############################################################################
# MODULE 1: Environment & Dependency Management
# ##############################################################################

init_setup <- function() {
  cat("\n>>> MODULE 1: Loading Packages\n")
  library(ggplot2)
  library(dplyr)
  library(tidyr)
}

load_data <- function() {
  cat(">>> MODULE 1: Loading Empirical Results\n")
  if (!file.exists("phase3a_empirical_expr.rds") | !file.exists("phase3a_empirical_meth.rds")) {
    stop("Run Scripts 01 and 02 first to generate the .rds files.")
  }
  return(list(
    expr = readRDS("phase3a_empirical_expr.rds"),
    meth = readRDS("phase3a_empirical_meth.rds")
  ))
}

# ##############################################################################
# MODULE 2: Data Preparation
# ##############################################################################

prepare_data <- function(data) {
  cat(">>> MODULE 2: Preparing Data Matrix\n")
  
  # Extract lists (FIX: Added Expression A)
  expr_A <- data$expr$lasso_compare$Strategy_A_Genes
  expr_B <- data$expr$lasso_compare$Strategy_B_Genes
  meth_A <- data$meth$results_A$gene_level$active_features
  meth_B <- data$meth$results_B$gene_level$active_features
  
  # Get all unique genes across all sets
  all_genes <- unique(c(expr_A, expr_B, meth_A, meth_B))
  
  # Create a binary matrix (1 if selected, 0 if not)
  df <- data.frame(Gene = all_genes)
  df$Expression_A   <- ifelse(df$Gene %in% expr_A, 1, 0)
  df$Expression_B   <- ifelse(df$Gene %in% expr_B, 1, 0)
  df$Methylation_A  <- ifelse(df$Gene %in% meth_A, 1, 0)
  df$Methylation_B  <- ifelse(df$Gene %in% meth_B, 1, 0)
  
  # Reshape to long format for ggplot
  df_long <- df %>%
    pivot_longer(cols = c("Expression_A", "Expression_B", "Methylation_A", "Methylation_B"),
                 names_to = "Strategy",
                 values_to = "Selected") %>%
    group_by(Gene) %>%
    mutate(Total_Selections = sum(Selected)) %>%
    ungroup() %>%
    mutate(
      # Order strategies logically
      Strategy = factor(Strategy, levels = c("Methylation_A", "Methylation_B", "Expression_A", "Expression_B"))
    ) %>%
    arrange(Total_Selections, Gene) %>%
    mutate(Gene = factor(Gene, levels = unique(Gene)))
  
  return(df_long)
}

# ##############################################################################
# MODULE 3: Plot Generation & Export
# ##############################################################################

generate_plot <- function(df_long) {
  cat(">>> MODULE 3: Generating Dot Plot\n")
  
  # Create the plot
  p <- ggplot(df_long, aes(x = Strategy, y = Gene)) +
    geom_tile(color = "grey90", fill = "white", linewidth = 0.5) +
    
    geom_point(data = subset(df_long, Selected == 1), 
               aes(color = Strategy), size = 3) +
    
    scale_color_manual(values = c(
      "Methylation_A" = "#D55E00", # Orange
      "Methylation_B" = "#0072B2", # Blue
      "Expression_A"  = "#CC79A7", # Purple (For Strategy A)
      "Expression_B"  = "#009E73"  # Green
    )) +
    
    scale_x_discrete(labels = c(
      "Methylation_A" = "Methylation\nStrategy A",
      "Methylation_B" = "Methylation\nStrategy B",
      "Expression_A"  = "Expression\nStrategy A",
      "Expression_B"  = "Expression\nStrategy B"
    )) +
    
    labs(
      title = "LASSO Feature Selection Overlap",
      subtitle = "Dots indicate the selection of a gene in the respective strategy",
      x = NULL,
      y = "Candidate Gene"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.y = element_text(face = "italic", size = 10),
      axis.text.x = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  # Note: Increased width to 8 to accommodate the 4th column
  ggsave("Figure3_Feature_Selection_DotPlot.pdf", p, width = 8, height = 8)
  cat("Saved Figure3_Feature_Selection_DotPlot.pdf successfully.\n")
}

# ##############################################################################
# MASTER EXECUTION BLOCK
# ##############################################################################

run_pipeline <- function() {
  cat("========================================================\n")
  cat("Starting Dot Plot Generation\n")
  cat("========================================================\n")
  
  init_setup()
  data <- load_data()
  df_long <- prepare_data(data)
  generate_plot(df_long)
  
  cat("\n========================================================\n")
  cat("Dot Plot Generation Complete!\n")
  cat("========================================================\n")
}

# Execute
run_pipeline()