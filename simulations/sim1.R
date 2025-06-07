# check, install, and load necessary packages
packages <- c("grf", "dplyr", "data.table", "tidyr")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))


# ===========================================================
# SIMULATION SCRIPT 1: Contributions to variance in Double ML
# -----------------------------------------------------------
# Runs a simulation to analyze variance sources in Double Machine Learning.
# Varies sample size, number of trees, and folds for cross-fitting.
# Decomposes total variance into sample, cross-fitting, and learner components.
# ===========================================================



# load functions
source("functions/generate_DGP_1.R")
source("functions/generate_DGP_2.R")
source("functions/get_estimates_sim1.R")
source("functions/variance_contributions.R")



set.seed(420)

##### Simulation Parameters #####
n_obs_values <- c(100, 200, 400, 800)
num_trees_values <- c(1, 10, 100, 500)
n_fold_values <- c(2, 5, 10)


# choose DGP function
# dgp_function <- generate_DGP_2
dgp_function <- generate_DGP_1

# used for output file name
dgp_name <- ifelse(identical(dgp_function, generate_DGP_1), "DGP1", "DGP2")

n_rep_dgp <- 100   # sampling repetitions
n_cv_rep <- 10     # cross-Fitting repetitions
n_rf_reps <- 10    # learner repetitions per cross-fitting round

# create parameter grid
param_grid <- expand.grid(
  n_obs = n_obs_values,
  n_trees = num_trees_values,
  n_folds = n_fold_values
)

final_results_all <- list()

##### Simulation Loop #####
for (i in 1:nrow(param_grid)) {
  n_obs <- param_grid$n_obs[i]
  n_trees <- param_grid$n_trees[i]
  n_folds <- param_grid$n_folds[i]
  
  message("Running: n_obs = ", n_obs, ", trees = ", n_trees, ", folds = ", n_folds)
  
  results_list_dml <- list()
  
  # sample repetitions
  for (dgp in 1:n_rep_dgp) {
    if (dgp %% 10 == 0) message("  DGP iteration: ", dgp)
    
    # generate data from chosen DGP
    data <- dgp_function(n_obs)
    Y <- data$y
    X <- data$x
    W <- data$w
    n <- length(Y)
    
    # cross-Fitting repetitions
    for (cv_rep in 1:n_cv_rep) {
      fold <- sample(rep(1:n_folds, length.out = n))
      exhat <- rep(NA, n)  
      mxhat <- rep(NA, n) 
      
      # learner repetitions
      for (rf_rep in 1:n_rf_reps) {
        for (i_fold in 1:n_folds) {
          
          # skip folds with only 1 treatment group
          if (length(unique(W[fold == i_fold])) < 2) next
          
          # fit models on training data
          rfe <- regression_forest(X[fold != i_fold, ], W[fold != i_fold], num.trees = n_trees)
          rfm <- regression_forest(X[fold != i_fold, ], Y[fold != i_fold], num.trees = n_trees)
          
          # predict on held-out fold
          exhat[fold == i_fold] <- predict(rfe, X[fold == i_fold, , drop = FALSE])$predictions
          mxhat[fold == i_fold] <- predict(rfm, X[fold == i_fold, , drop = FALSE])$predictions
        }
        
        valid_idx <- which(!is.na(exhat))
        if (length(valid_idx) == 0) next
        
        # truncate extreme propensity scores
        exhat[valid_idx] <- pmax(pmin(exhat[valid_idx], 1 - 1e-2), 1e-2)
        
        # compute orthogonalized score components
        psi_a <- -(W[valid_idx] - exhat[valid_idx])^2
        psi_b <- (Y[valid_idx] - mxhat[valid_idx]) * (W[valid_idx] - exhat[valid_idx])
        
        # Double ML estimation
        dml_result <- dml_inference_sim1(psi_a, psi_b)
        
        # store results
        results_list_dml[[length(results_list_dml) + 1]] <- data.table(
          dgp = dgp,
          cv_rep = cv_rep,
          rf_rep = rf_rep,
          theta_hat = dml_result[1],
          se_hat = dml_result[2],
          t_value = dml_result[3],
          p_value = dml_result[4]
        )
      }
    }
  }
  
  # get contributions of each layer
  final_results_dml <- rbindlist(results_list_dml)
  variance_results <- variance_decomposition_sim1(final_results_dml)
  
  # save combination of parameters and results
  final_results_all[[length(final_results_all) + 1]] <- cbind(
    data.table(n_obs = n_obs, n_trees = n_trees, n_folds = n_folds),
    variance_results
  )
}

# combine all results into one table
final_variance_results <- rbindlist(final_results_all)

# save with dynamic filename
fwrite(final_variance_results, file = paste0("sim1_results_", dgp_name, ".csv"))

