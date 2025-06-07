# check, install, and load necessary packages
packages <- c("xgboost", "caret")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))

#' Hyperparameter Tuning for XGBoost using 5-Fold Cross-Validation
#'
#' Performs grid search over a set of XGBoost hyperparameter combinations using 5-fold cross-validation.
#' Includes early stopping and tracks the average number of boosting iterations across folds.
#'
#' @param x_train A data.frame or matrix containing the training features.
#' @param y_train A numeric vector containing the training outcome.
#' @param param_grid A data.frame containing the hyperparameter combinations to evaluate.
#'   expected columns: eta, max_depth, nrounds.
#' @param metric A string specifying the evaluation metric. Either "mse" (mean squared error) or "log_loss" (logarithmic loss).
#'
#' @return A named list containing the best-performing hyperparameter values, including the average optimal number of boosting rounds.
tune_xgb <- function(x_train, y_train, param_grid, metric = "mse") {
  best_params <- NULL
  best_score <- Inf  
  best_nrounds <- NULL  
  
  # create 5 CV folds
  folds_tune <- createFolds(y_train, k = 5, list = TRUE)
  
  # loop through all combinations in the grid
  for (params in seq_len(nrow(param_grid))) {
    # extract current parameter combination
    eta <- param_grid[params, "eta"]
    max_depth <- param_grid[params, "max_depth"]
    nrounds <- param_grid[params, "nrounds"]
    
    fold_scores <- c()     
    fold_nrounds <- c()    
    
    # loop over folds
    for (fold_tune in folds_tune) {
      # split data into training and validation
      x_train_fold <- x_train[-fold_tune, ]
      y_train_fold <- y_train[-fold_tune]
      x_valid_fold <- x_train[fold_tune, ]
      y_valid_fold <- y_train[fold_tune]
      
      # set objective and eval metric
      if (metric == "mse") {
        objective <- "reg:squarederror"
        eval_metric <- "rmse"
      } else if (metric == "log_loss") {
        objective <- "binary:logistic"
        eval_metric <- "logloss"
      }
      
      # prepare xgboost parameter list
      params_list <- list(
        objective = objective,
        eta = eta,
        max_depth = max_depth,
        lambda = 1,
        alpha = 0.1,
        eval_metric = eval_metric
      )
      
      # convert to DMatrix
      dtrain <- xgb.DMatrix(data = as.matrix(x_train_fold), label = y_train_fold)
      dvalid <- xgb.DMatrix(data = as.matrix(x_valid_fold), label = y_valid_fold)
      
      # train xgboost model with early stopping
      xgb_model <- xgb.train(
        params = params_list,
        data = dtrain,
        nrounds = nrounds,
        watchlist = list(train = dtrain, valid = dvalid),
        early_stopping_rounds = 8,
        verbose = 0
      )
      
      # predict on validation set
      preds <- predict(xgb_model, newdata = dvalid, iteration_range = xgb_model$best_iteration)
      
      # evaluate performance
      if (metric == "mse") {
        score <- mean((y_valid_fold - preds)^2)
      } else if (metric == "log_loss") {
        preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)  # avoid log(0)
        score <- -mean(y_valid_fold * log(preds) + (1 - y_valid_fold) * log(1 - preds))
      }
      
      fold_scores <- c(fold_scores, score)
      fold_nrounds <- c(fold_nrounds, xgb_model$best_iteration)
    }
    
    # compute average score and rounds over folds
    avg_score <- mean(fold_scores)
    avg_nrounds <- floor(mean(fold_nrounds))
    
    # update best parameters if improved
    if (avg_score < best_score) {
      best_score <- avg_score
      best_params <- list(
        eta = eta,
        max_depth = max_depth,
        nrounds = avg_nrounds
      )
    }
  }
  
  return(best_params)
}
