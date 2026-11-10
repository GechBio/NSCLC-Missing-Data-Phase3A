# ==============================================================================
# Script 04b: Generate Ground Truth Figure (Kaplan-Meier) with 95% CI
#
# Author:   Getachew Derbew Belay
# Date:     November 2026
# Platform: R 4.5.2
#
# Purpose:  Generate the Kaplan-Meier survival curves for the calibrated 
#           simulation ground truth. Automatically extracts the Hazard Ratios 
#           and 95% CIs and annotates them directly onto the figure.
#
# Input:    phase3a_sim_ground_truth.rds
# Output:   Figure5_Simulation_Ground_Truth_KM.pdf
# ==============================================================================

# ##############################################################################
# MODULE 1: Setup & Data Loading
# ##############################################################################

cat("\n>>> MODULE 1: Loading Data and Packages\n")

library(ggplot2)
library(survival)
library(survminer)

if (!file.exists("phase3a_sim_ground_truth.rds")) {
  stop("Run Script 04 first to generate the ground truth.")
}

gt_data <- readRDS("phase3a_sim_ground_truth.rds")
df <- gt_data$ground_truth_data

# Ensure stage_group is a factor with the correct order
df$stage_group <- factor(df$stage_group, levels = c("Early", "Advanced", "Unknown"))

# ##############################################################################
# MODULE 2: Calculate HRs and 95% CIs
# ##############################################################################

cat(">>> MODULE 2: Calculating Hazard Ratios and 95% CIs\n")

# Fit the Cox model to extract the exact HRs and 95% CIs
cox_fit <- coxph(Surv(sim_time, sim_event) ~ stage_group, data = df)

# Extract Hazard Ratios
hr_adv <- exp(coef(cox_fit))["stage_groupAdvanced"]
hr_unk <- exp(coef(cox_fit))["stage_groupUnknown"]

# Extract 95% Confidence Intervals
ci <- exp(confint(cox_fit))
ci_adv_lower <- ci["stage_groupAdvanced", 1]
ci_adv_upper <- ci["stage_groupAdvanced", 2]
ci_unk_lower <- ci["stage_groupUnknown", 1]
ci_unk_upper <- ci["stage_groupUnknown", 2]

# Create the text annotations for the plot
text_adv <- sprintf("Advanced HR: %.2f (95%% CI: %.2f-%.2f)", hr_adv, ci_adv_lower, ci_adv_upper)
text_unk <- sprintf("Unknown HR: %.2f (95%% CI: %.2f-%.2f)", hr_unk, ci_unk_lower, ci_unk_upper)
annotation_text <- paste("vs. Early Stage:\n", text_adv, "\n", text_unk)

# ##############################################################################
# MODULE 3: Generate Kaplan-Meier Plot
# ##############################################################################

cat(">>> MODULE 3: Generating Kaplan-Meier Plot\n")

# Create the Kaplan-Meier fit
km_fit <- survfit(Surv(sim_time, sim_event) ~ stage_group, data = df)

# Generate the plot
fig5 <- ggsurvplot(
  km_fit,
  data = df,
  size = 1.2,
  palette = c("#2E9FDF", "#E7B800", "#FC4E07"), # Blue, Gold, Red
  conf.int = FALSE,
  pval = TRUE,
  risk.table = FALSE, 
  xlab = "Time (Days)",
  ylab = "Overall Survival Probability",
  legend.title = "True Stage",
  legend.labs = c("Early", "Advanced", "Unknown (Latent)"),
  ggtheme = theme_classic(),
  title = "Simulated Ground Truth Survival Curves (N=867)"
)

# Add the HR and 95% CI annotation to the plot
fig5$plot <- fig5$plot + 
  annotate("text", x = 3000, y = 0.85, 
           label = annotation_text, 
           hjust = 0, size = 4, fontface = "bold", color = "black")

# ##############################################################################
# MODULE 4: Save Figure
# ##############################################################################

cat(">>> MODULE 4: Saving Figure\n")

ggsave(
  filename = "Figure5_Simulation_Ground_Truth_KM.pdf", 
  plot = fig5$plot, 
  width = 8, 
  height = 6,
  device = "pdf"
)

cat("Saved Figure5_Simulation_Ground_Truth_KM.pdf successfully with 95% CI annotations.\n")