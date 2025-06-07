# check, install, and load necessary packages
packages <- c("grf", "xgboost", "glmnet", "data.table", "caret")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))


# ===========================================================
# SIMULATION SCRIPT 2 — Repeated cross-fitting on a single sample
# -----------------------------------------------------------
# For one synthetic dataset: repeats cross-fitting with multiple learners,
# tunes hyperparameters, collects all model predictions and saves results.
# Design: one sample from dgp, n = 400, 500 cross-fitting reps (5 folds)
# focus: aggregation/learner/estimator variation on identical data
# ===========================================================


# load functions
source("functions/generate_DGP_1.R")
source("functions/generate_DGP_2.R")
source("functions/step_function.R")
source("functions/tune_lasso.R")
source("functions/tune_rf.R")
source("functions/tune_xgb.R")
source("functions/transform_x_lasso.R")
source("functions/learner_methods.R")
source("functions/hyperparams_grid.R")






# config
set.seed(420)
n_cv_rep <- 500  
nfolds <- 5      
n_obs <- 400     
learners <- c("xgboost", "random_forest", "lasso") 
dgp_choice <- "DGP_1"
# dgp_choice <- "DGP_2" 


# initialize result container
all_results_list_all <- list()

# generate one sample
data <- switch(dgp_choice,
               "DGP_1" = generate_DGP_1(n_obs),
               "DGP_2" = generate_DGP_2(n_obs))

y <- data$y
x <- data$x
w <- data$w

# cross-validation fold indices for each repetition
folds_list <- replicate(n_cv_rep, createFolds(1:n_obs, k = nfolds, list = FALSE), simplify = FALSE)

# simulation loop
for (learner in learners) {
  cat("Current learner:", learner, "\n")
  
  methods <- learner_methods[[learner]]
  x_input <- if (learner == "lasso") transform_x_lasso(x) else x
  
  # hyperparameter tuning per learner
  best_params_w   <- methods$tune(x_input, w, ifelse(learner == "lasso", "binomial", "log_loss"))
  best_params_y   <- methods$tune(x_input, y, ifelse(learner == "lasso", "gaussian", "mse"))
  best_params_y_w0 <- methods$tune(x_input[w == 0, ], y[w == 0], ifelse(learner == "lasso", "gaussian", "mse"))
  best_params_y_w1 <- methods$tune(x_input[w == 1, ], y[w == 1], ifelse(learner == "lasso", "gaussian", "mse"))
  best_params_str <- function(p) paste(unlist(p), collapse = "_")
  
  # initialize psi arrays
  exhat_mat   <- matrix(NA, nrow = n_obs, ncol = n_cv_rep)
  mxhat_mat   <- matrix(NA, nrow = n_obs, ncol = n_cv_rep)
  mwhat0_mat  <- matrix(NA, nrow = n_obs, ncol = n_cv_rep)
  mwhat1_mat  <- matrix(NA, nrow = n_obs, ncol = n_cv_rep)
  
  for (cv_rep in 1:n_cv_rep) {
    if (cv_rep %% 50 == 0) cat("Cross-Fitting Repetition:", cv_rep, "\n")
    
    fold <- folds_list[[cv_rep]]
    exhat <- mxhat <- mwhat0 <- mwhat1 <- rep(NA, n_obs)
    
    for (i in 1:nfolds) {
      idx_train <- which(fold != i)
      idx_test  <- which(fold == i)
      
      # model training and prediction
      if (learner == "xgboost") {
        model_e  <- methods$train(x_input[idx_train, ], w[idx_train], best_params_w, "binary:logistic")
        model_m  <- methods$train(x_input[idx_train, ], y[idx_train], best_params_y, "reg:squarederror")
        model_m0 <- methods$train(x_input[idx_train & w[idx_train]==0, ], y[idx_train & w[idx_train]==0], best_params_y_w0, "reg:squarederror")
        model_m1 <- methods$train(x_input[idx_train & w[idx_train]==1, ], y[idx_train & w[idx_train]==1], best_params_y_w1, "reg:squarederror")
        
        exhat[idx_test]  <- methods$predict(model_e,  x_input[idx_test, ])
        mxhat[idx_test]  <- methods$predict(model_m,  x_input[idx_test, ])
        mwhat0[idx_test] <- methods$predict(model_m0, x_input[idx_test, ])
        mwhat1[idx_test] <- methods$predict(model_m1, x_input[idx_test, ])
        
      } else if (learner == "lasso") {
        model_e  <- methods$train(x_input[idx_train, ], w[idx_train], best_params_w, "binomial")
        model_m  <- methods$train(x_input[idx_train, ], y[idx_train], best_params_y, "gaussian")
        model_m0 <- methods$train(x_input[idx_train & w[idx_train]==0, ], y[idx_train & w[idx_train]==0], best_params_y_w0, "gaussian")
        model_m1 <- methods$train(x_input[idx_train & w[idx_train]==1, ], y[idx_train & w[idx_train]==1], best_params_y_w1, "gaussian")
        
        exhat[idx_test]  <- methods$predict(model_e,  x_input[idx_test, ], best_params_w,  "response")
        mxhat[idx_test]  <- methods$predict(model_m,  x_input[idx_test, ], best_params_y,  "link")
        mwhat0[idx_test] <- methods$predict(model_m0, x_input[idx_test, ], best_params_y_w0, "link")
        mwhat1[idx_test] <- methods$predict(model_m1, x_input[idx_test, ], best_params_y_w1, "link")
        
      } else if (learner == "random_forest") {
        model_e  <- methods$train(x_input[idx_train, ], w[idx_train], best_params_w)
        model_m  <- methods$train(x_input[idx_train, ], y[idx_train], best_params_y)
        model_m0 <- methods$train(x_input[idx_train & w[idx_train]==0, ], y[idx_train & w[idx_train]==0], best_params_y_w0)
        model_m1 <- methods$train(x_input[idx_train & w[idx_train]==1, ], y[idx_train & w[idx_train]==1], best_params_y_w1)
        
        exhat[idx_test]  <- methods$predict(model_e,  x_input[idx_test, ])
        mxhat[idx_test]  <- methods$predict(model_m,  x_input[idx_test, ])
        mwhat0[idx_test] <- methods$predict(model_m0, x_input[idx_test, ])
        mwhat1[idx_test] <- methods$predict(model_m1, x_input[idx_test, ])
      }
    }
    
    # store predictions for cross-fitting repetition
    exhat_mat[, cv_rep]   <- exhat
    mxhat_mat[, cv_rep]   <- mxhat
    mwhat0_mat[, cv_rep]  <- mwhat0
    mwhat1_mat[, cv_rep]  <- mwhat1
  }
  
  # store results
  all_results_list_all[[length(all_results_list_all) + 1]] <- data.table(
    learner = learner,
    best_params_w = best_params_str(best_params_w),
    best_params_y = best_params_str(best_params_y),
    best_params_y_w0 = best_params_str(best_params_y_w0),
    best_params_y_w1 = best_params_str(best_params_y_w1),
    x = list(x),
    y = list(y),
    w = list(w),
    exhat = list(exhat_mat),
    mxhat = list(mxhat_mat),
    mwhat0 = list(mwhat0_mat),
    mwhat1 = list(mwhat1_mat)
  )
}



# combine and save results
final_results <- rbindlist(all_results_list_all)
output_filename <- paste0("sim2_results_", dgp_choice, ".rds")
saveRDS(final_results, output_filename, compress = "xz")

