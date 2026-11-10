# Phase 3A: Missing Data Strategies in NSCLC

This repository contains the complete R pipeline to reproduce the findings in the manuscript: *"Evaluating Missing-Data Strategies for Pathologic Stage in High-Dimensional NSCLC Prognostic Modeling: A Calibrated Simulation and Empirical Study."*

## Contents
* `00_base_setup_phase3a.R` to `08_master_pipeline.R`: The modular R scripts.
* `renv.lock`: The frozen R environment containing exact package versions.
* Generated figures (`.pdf`) and tables (`.csv`).

## How to Reproduce the Results
1. Ensure R (v4.5.2) is installed.
2. Restore the package environment:
   ```R
   install.packages("renv")
   renv::restore()