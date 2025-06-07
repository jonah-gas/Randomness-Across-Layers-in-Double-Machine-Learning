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
# SIMULATION SCRIPT 3: Inferrnce over repeated samples in Double ML 
# -----------------------------------------------------------
# Repeatedly draws samples from dgp, runs cross-fitting and estimation
# with multiple learners, collects psi for PLR and AIPW, and saves results.
# Design: 500 samples per dgp, n = 400, 100 cross-fitting reps (5 folds).
# Focus: disentangle systematic patterns from random noise in estimation
# ===========================================================


# load functions
source("functions/generate_DGP_1.R")
source("functions/generate_DGP_2.R")
source("functions/step_function.R")
source("functions/tune_lasso.R")
source("functions/tune_rf.R")
source("functions/tune_xgb.R")
source("functions/transform_x_lasso.R")
source("functions/hyperparams_grid.R")
source("functions/learner_methods.R")
source("functions/get_estimates.R")

# config
set.seed(420)
n_cv_rep <- 100
nfolds <- 5
n_obs <- 400
n_sim <- 500
trunc <- 0.01
learners <- c("lasso", "random_forest", "xgboost")
dgp_choice <- "DGP_1"
#dgp_choice <- "DGP_2"
n_learners <- length(learners)

# initialize result container
estimates_list <- list()

# simulation loop
for (sim in 1:n_sim) {
  cat("\n===== simulation", sim, "=====\n")
  t_start <- Sys.time()
  
  # generate data
  data <- switch(dgp_choice,
                 "DGP_1" = generate_DGP_1(n_obs),
                 "DGP_2" = generate_DGP_2(n_obs))
  y <- data$y
  x <- data$x
  w <- data$w
  folds_list <- replicate(n_cv_rep, createFolds(1:n_obs, k = nfolds, list = FALSE), simplify = FALSE)
  
  # temporary storage per simulation
  sim_result <- list(simulation = sim, learners = list())
  
  for (l in seq_along(learners)) {
    learner <- learners[l]
    cat(" → learner:", learner, "\n")
    methods <- learner_methods[[learner]]
    x_input <- if (learner == "lasso") transform_x_lasso(x) else x
    
    # hyperparameter tuning per learner and simulation
    best_params_w    <- methods$tune(x_input, w, ifelse(learner == "lasso", "binomial", "log_loss"))
    best_params_y    <- methods$tune(x_input, y, ifelse(learner == "lasso", "gaussian", "mse"))
    best_params_y_w0 <- methods$tune(x_input[w == 0, ], y[w == 0], ifelse(learner == "lasso", "gaussian", "mse"))
    best_params_y_w1 <- methods$tune(x_input[w == 1, ], y[w == 1], ifelse(learner == "lasso", "gaussian", "mse"))
    
    # initialize psi arrays
    psi_a_plr_tmp <- matrix(NA, nrow = n_obs, ncol = n_cv_rep)
    psi_b_plr_tmp <- matrix(NA, nrow = n_obs, ncol = n_cv_rep)
    psi_a_irm_tmp <- matrix(NA, nrow = n_obs, ncol = n_cv_rep)
    psi_b_irm_tmp <- matrix(NA, nrow = n_obs, ncol = n_cv_rep)
    
    
    for (cv_rep in 1:n_cv_rep) {
      if (cv_rep %% 50 == 0) cat("   repetition", cv_rep, "\n")
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
      
      # compute influence functions
      keep <- (exhat > trunc) & (exhat < 1 - trunc)
      psi_a_plr_tmp[keep, cv_rep] <- -(w[keep] - exhat[keep])^2
      psi_b_plr_tmp[keep, cv_rep] <- (y[keep] - mxhat[keep]) * (w[keep] - exhat[keep])
      psi_a_irm_tmp[keep, cv_rep] <- -1
      psi_b_irm_tmp[keep, cv_rep] <- mwhat1[keep] - mwhat0[keep] +
        w[keep] * (y[keep] - mwhat1[keep]) / exhat[keep] -
        (1 - w[keep]) * (y[keep] - mwhat0[keep]) / (1 - exhat[keep])
    }
    
    # collect results
    est_plr <- get_estimates(array(psi_a_plr_tmp, dim = c(n_obs, n_cv_rep, 1)),
                             array(psi_b_plr_tmp, dim = c(n_obs, n_cv_rep, 1)))
    est_irm <- get_estimates(array(psi_a_irm_tmp, dim = c(n_obs, n_cv_rep, 1)),
                             array(psi_b_irm_tmp, dim = c(n_obs, n_cv_rep, 1)))
    
    # save in list
    estimates_list[[length(estimates_list) + 1]] <- data.frame(
      sim = sim,
      learner = learner,
      method = "PLR",
      as.data.frame(est_plr)
    )
    estimates_list[[length(estimates_list) + 1]] <- data.frame(
      sim = sim,
      learner = learner,
      method = "IRM",
      as.data.frame(est_irm)
    )
  }
  
  t_end <- Sys.time()
  cat("⏱ duration:", round(difftime(t_end, t_start, units = "mins"), 2), "min\n")
}

# save full result
estimates_df <- do.call(rbind, estimates_list)
saveRDS(estimates_df, file = paste0("sim3_results_", dgp_choice, ".rds"))