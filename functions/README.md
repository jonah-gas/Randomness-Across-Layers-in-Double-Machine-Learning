# functions/

This folder contains all R functions and helper scripts used in the project.

**Overview of contents:**
- **Data Generating Processes:**  
  Functions for creating synthetic data sets (e.g., `generate_DGP_1.R`, `generate_DGP_2.R`)
- **Feature Engineering:**  
  Helper functions for feature expansion and transformation (`transform_x_lasso.R`, `step_function.R`)
- **Model Training and Tuning:**  
  Functions and lists for model hyperparameter tuning and training (`tune_lasso.R`, `tune_rf.R`, `tune_xgb.R`, `learner_methods.R`, `hyperparams_grid.R`)
- **Estimation & Inference:**  
  Functions for treatment effect estimation and variance decomposition (`get_estimates.R`, `variance_contributions.R`, `get_estimates_sim1.R`)

All scripts automatically check and load their required R packages.
