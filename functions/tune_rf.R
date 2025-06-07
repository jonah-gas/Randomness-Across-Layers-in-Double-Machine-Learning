# check, install, and load necessary packages
packages <- c("grf", "caret")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))

#' Hyperparameter Tuning for Random Forest using 5-Fold Cross-Validation
#'
#' Performs grid search over a set of random forest hyperparameters using 5-fold cross-validation. 
#' For each parameter combination, cross-validated performance is evaluated using the specified error metric.
#' The function returns the hyperparameter set that minimizes the average cross-validation error.
#'
#' @param x_train A data.frame or matrix containing the training features.
#' @param y_train A numeric vector containing the training outcome.
#' @param param_grid A data.frame containing the hyperparameter combinations to evaluate.
#'   expected columns: num.trees, sample.fraction, min.node.size, mtry.
#' @param metric A string specifying the evaluation metric. Either "mse" (mean squared error) or "log_loss" (logarithmic loss).
#'
#' @return A named list containing the best-performing hyperparameter values.
#' 
tune_rf <- function(x_train, y_train, param_grid, metric = "mse") {
  best_params <- NULL
  best_score <- Inf  
  
  # create 5-fold cross-validation splits
  folds <- createFolds(y_train, k = 5, list = TRUE)
  
  # loop over each hyperparameter combination in the grid
  for (params in seq_len(nrow(param_grid))) {
    # extract current set of hyperparameters
    num.trees <- param_grid[params, "num.trees"]
    sample.fraction <- param_grid[params, "sample.fraction"]
    min.node.size <- param_grid[params, "min.node.size"]
    mtry <- param_grid[params, "mtry"]
    
    fold_scores <- c()  
    
    # cross-validation loop
    for (fold in folds) {
      # define training and validation data for the current fold
      x_train_fold <- as.matrix(x_train[-fold, ])
      y_train_fold <- y_train[-fold]
      x_valid_fold <- as.matrix(x_train[fold, ])
      y_valid_fold <- y_train[fold]
      
      # train random forest on training fold
      rf <- regression_forest(
        x_train_fold, y_train_fold,
        num.trees = num.trees,
        sample.fraction = sample.fraction,
        min.node.size = min.node.size,
        mtry = mtry,
        honesty = FALSE,
        ci.group.size = 1
      )
      
      # predict on validation fold
      preds <- predict(rf, x_valid_fold)$predictions
      
      # compute performance metric
      if (metric == "mse") {
        score <- mean((y_valid_fold - preds)^2)  # mean squared error
        fold_scores <- c(fold_scores, score)
      } 
      
      if (metric == "log_loss") {
        # ensure predictions are within (0,1) to avoid log(0)
        preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)
        score <- -mean(y_valid_fold * log(preds) + (1 - y_valid_fold) * log(1 - preds))  # logarithmic loss
        fold_scores <- c(fold_scores, score)
      }
    }
    
    # compute average score across folds
    avg_score <- mean(fold_scores)
    
    # update best parameters if current configuration improves the score
    if (avg_score < best_score) {
      best_score <- avg_score
      best_params <- list(
        num.trees = num.trees,
        sample.fraction = sample.fraction,
        min.node.size = min.node.size,
        mtry = mtry
      )
    }
  }
  
  return(best_params)
}
