# Phase 3A: Missing Data Strategies in NSCLC

This repository contains the complete R pipeline to reproduce the findings in the manuscript: *"Evaluating Missing-Data Strategies for Pathologic Stage in High-Dimensional NSCLC Prognostic Modeling: A Calibrated Simulation and Empirical Study."*

## Contents

* `00_base_setup_phase3a.R` to `08_master_pipeline.R`: The modular R scripts (pipeline steps).
* `09_clinical_cohort_disambiguation.R`: Supplementary analysis (Tables S1 and S2 — KM median vs observation time; internal consistency checks).
* `10_ipw_controlled_test.R`: Supplementary analysis (Table S3 — IPW pure-MAR control).
* `simulation_helpers.R`: Helper library (parameterized ground-truth functions) sourced by Script 10. Not a pipeline step.
* `renv.lock`: The frozen R environment containing exact package versions.
* Generated figures (`.pdf`) and tables (`.csv`).

## Environment Setup

1. Ensure R (v4.5.2) is installed.

2. Restore the package environment:
   ```R
   install.packages("renv")
   renv::restore()