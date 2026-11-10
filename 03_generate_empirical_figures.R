# ==============================================================================
# Script 03: Generate Empirical Arm Figures (C-index & HR)
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
#
# Purpose:  Load the results from Scripts 01 and 02, and generate required figures:
#           Figure 2: C-index Dumbbell Plot (Expression Models)
#           Figure 4: Annotated Forest Plot (HR per IQR for Methylation)
#
# Input:    phase3a_empirical_expr.rds, phase3a_empirical_meth.rds
# Output:   Figure2_Cindex_Comparison.pdf
#           Figure4_Methylation_HR_Annotated.pdf
# ==============================================================================

# ##############################################################################
# MODULE 1: Setup & Data Loading
# ##############################################################################

cat("\n>>> MODULE 1: Loading Data and Packages\n")

library(ggplot2)
library(dplyr)
library(tidyr)

if (!file.exists("phase3a_empirical_expr.rds") | !file.exists("phase3a_empirical_meth.rds")) {
  stop("Run Scripts 01 and 02 first to generate the .rds files.")
}

expr_data <- readRDS("phase3a_empirical_expr.rds")
meth_data <- readRDS("phase3a_empirical_meth.rds")


# ##############################################################################
# MODULE 2: Figure 2 - C-index Dumbbell Plot (Expression)
# ##############################################################################

cat("\n>>> MODULE 2: Generating Figure 2 (C-index Comparison)...\n")

cindex_df <- data.frame(
  Model = c("Clinical", "Stepwise Cox", "RSF", "Elastic Net", "LASSO", "Ridge"),
  Strategy_A = expr_data$cindex_table$Strategy_A_Cindex,
  Strategy_B = expr_data$cindex_table$Strategy_B_Cindex
)

# Reshape to long format for ggplot
cindex_long <- cindex_df %>%
  pivot_longer(cols = c("Strategy_A", "Strategy_B"), 
               names_to = "Strategy", values_to = "Cindex") %>%
  mutate(Strategy = ifelse(Strategy == "Strategy_A", "Strategy A (Complete-Case)", "Strategy B (Missing-Retained)"))

fig2 <- ggplot(cindex_long, aes(x = Cindex, y = reorder(Model, Cindex))) +
  # Line connecting the two points
  geom_line(aes(group = Model), color = "gray60", linewidth = 1) +
  # Points for Strategy A and B
  geom_point(aes(color = Strategy), size = 4) +
  # Random baseline line
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.6) +
  scale_color_manual(values = c("Strategy A (Complete-Case)" = "#D55E00", 
                                "Strategy B (Missing-Retained)" = "#0072B2")) +
  labs(
    title = "Validation C-index: Strategy A vs. Strategy B",
    subtitle = "Expression Cohort (N=744 vs. N=867)",
    x = "Validation C-index",
    y = NULL,
    color = "Strategy"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 11)
  )

ggsave("Figure2_Cindex_Comparison.pdf", fig2, width = 7, height = 5)
cat("Saved Figure2_Cindex_Comparison.pdf\n")


# ##############################################################################
# MODULE 3: Figure 4 - Annotated Forest Plot (Methylation HR)
# ##############################################################################

cat("\n>>> MODULE 3: Generating Figure 4 (Annotated HR Forest Plot)...\n")

hr_df <- data.frame(
  Strategy = c("Strategy A\n(Complete-Case)", "Strategy B\n(Missing-Retained)"),
  HR = c(meth_data$results_A$gene_level$hr_iqr, meth_data$results_B$gene_level$hr_iqr),
  Lower = c(meth_data$results_A$gene_level$hr_ci[1], meth_data$results_B$gene_level$hr_ci[1]),
  Upper = c(meth_data$results_A$gene_level$hr_ci[2], meth_data$results_B$gene_level$hr_ci[2]),
  Pval = c(meth_data$results_A$gene_level$p_value, meth_data$results_B$gene_level$p_value)
)

# Create the text labels for the plot
hr_df$Label <- sprintf("HR: %.2f (95%% CI: %.2f-%.2f)\np = %.3f", 
                       hr_df$HR, hr_df$Lower, hr_df$Upper, hr_df$Pval)

fig4 <- ggplot(hr_df, aes(x = HR, y = Strategy, color = Strategy)) +
  # Add vertical reference line at HR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  
  # Add the point estimate and error bars
  geom_pointrange(aes(xmin = Lower, xmax = Upper), size = 1.2, linewidth = 1.2) +
  
  # Add the explicit text labels (HR, CI, and p-value)
  geom_text(aes(label = Label), 
            vjust = -1.2, # Position text above the points
            size = 4.5, 
            fontface = "bold",
            color = "black") +
  
  # Define colors
  scale_color_manual(values = c("Strategy A\n(Complete-Case)" = "#D55E00", 
                                "Strategy B\n(Missing-Retained)" = "#0072B2")) +
  
  # Adjust X-axis to give room for the text
  scale_x_continuous(limits = c(0.8, 1.9), breaks = seq(0.8, 1.9, by = 0.2)) +
  
  # Labels and Titles
  labs(
    title = "Methylation Risk Score: Adjusted HR per IQR",
    subtitle = "Adjusted for Age, Gender, Histology, and Stage (Full Cohort N=701)",
    x = "Hazard Ratio (95% Confidence Interval)",
    y = NULL
  ) +
  
  # Clean, publication-ready theme
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 15)),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10))
  )

ggsave("Figure4_Methylation_HR_Annotated.pdf", fig4, width = 8, height = 5)
cat("Saved Figure4_Methylation_HR_Annotated.pdf\n")

cat("\n========================================================\n")
cat("All empirical figures successfully generated as PDFs.\n")
cat("========================================================\n")