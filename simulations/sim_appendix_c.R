# check, install, and load necessary packages
packages <- c("grf", "data.table", "caret")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))



# ===========================================================
# SIMULATION SCRIPT Appendix: Inference over repeated samples 
# in Double ML for larger n
# -----------------------------------------------------------
# Additional Analysis similar to Simulation 3 for n = 2000,
# but only for DGP2, random forest, and AIPW. We reduce the
# number of samples to 100 and the cross-fitting repetitions
# to 10
# ===========================================================


# load functions
source("functions/generate_DGP_2.R")
source("functions/step_function.R")
source("functions/tune_rf.R")
source("functions/hyperparams_grid.R")
source("functions/learner_methods.R")
source("functions/get_estimates.R")

# config
set.seed(420)
n_cv_rep <- 10
nfolds <- 5
n_obs <- 2000
n_sim <- 100
trunc <- 0.01
learners <- c("random_forest")
dgp_choice <- "DGP_2"
n_learners <- length(learners)



# initialize result container
estimates_list <- list()


# simulation loop
for (sim in 1:n_sim) {
  cat("\n===== Simulation", sim, "=====\n")
  t_start <- Sys.time()
  
  # Daten generieren
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
    cat(" → Learner:", learner, "\n")
    methods <- learner_methods[[learner]]
    x_input <- if (learner == "lasso") transform_x_lasso(x) else x
    
    # hyperparameter tuning per learner and simulation
    best_params_w    <- methods$tune(x_input, w, ifelse(learner == "lasso", "binomial", "log_loss"))
    best_params_y_w0 <- methods$tune(x_input[w == 0, ], y[w == 0], ifelse(learner == "lasso", "gaussian", "mse"))
    best_params_y_w1 <- methods$tune(x_input[w == 1, ], y[w == 1], ifelse(learner == "lasso", "gaussian", "mse"))
    
    # initialize psi arrays
    psi_a_irm_tmp <- matrix(NA, nrow = n_obs, ncol = n_cv_rep)
    psi_b_irm_tmp <- matrix(NA, nrow = n_obs, ncol = n_cv_rep)
    
    
    for (cv_rep in 1:n_cv_rep) {
      if (cv_rep %% 50 == 0) cat("   Repetition", cv_rep, "\n")
      fold <- folds_list[[cv_rep]]
      exhat <- mxhat <- mwhat0 <- mwhat1 <- rep(NA, n_obs)
      
      for (i in 1:nfolds) {
        idx_train <- which(fold != i)
        idx_test  <- which(fold == i)
        # model training and prediction
        if (learner == "random_forest") {
          model_e  <- methods$train(x_input[idx_train, ], w[idx_train], best_params_w)
          model_m0 <- methods$train(x_input[idx_train & w[idx_train]==0, ], y[idx_train & w[idx_train]==0], best_params_y_w0)
          model_m1 <- methods$train(x_input[idx_train & w[idx_train]==1, ], y[idx_train & w[idx_train]==1], best_params_y_w1)
          
          exhat[idx_test]  <- methods$predict(model_e,  x_input[idx_test, ])
          mwhat0[idx_test] <- methods$predict(model_m0, x_input[idx_test, ])
          mwhat1[idx_test] <- methods$predict(model_m1, x_input[idx_test, ])
        }
      }
      
      # compute influence functions
      keep <- (exhat > trunc) & (exhat < 1 - trunc)
      psi_a_irm_tmp[keep, cv_rep] <- -1
      psi_b_irm_tmp[keep, cv_rep] <- mwhat1[keep] - mwhat0[keep] +
        w[keep] * (y[keep] - mwhat1[keep]) / exhat[keep] -
        (1 - w[keep]) * (y[keep] - mwhat0[keep]) / (1 - exhat[keep])
    }
    
    # collect results
    est_irm <- get_estimates(array(psi_a_irm_tmp, dim = c(n_obs, n_cv_rep, 1)),
                             array(psi_b_irm_tmp, dim = c(n_obs, n_cv_rep, 1)))
    
    estimates_list[[length(estimates_list) + 1]] <- data.frame(
      sim = sim,
      learner = learner,
      method = "IRM",
      as.data.frame(est_irm)
    )
  }
  
  
  t_end <- Sys.time()
  cat("⏱ Duration:", round(difftime(t_end, t_start, units = "mins"), 2), "min\n")
}

# --- Save full result ---

saveRDS(estimates_list, file = "sim_appendix_c_results.rds")

