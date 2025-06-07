# check, install, and load necessary packages
packages <- c("glmnet")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))

#' Hyperparameter Tuning for Lasso using 5-Fold Cross-Validation
#'
#' Performs 5-fold cross-validation to select the optimal penalty parameter (lambda)
#' for lasso regression using the glmnet package. Optionally includes polynomial and interaction features.
#'
#' @param x_train A data.frame or matrix containing the training features.
#' @param y_train A numeric vector containing the training outcome.
#' @param family A character string specifying the model family to be used in glmnet (e.g., "gaussian" or "binomial").
#' @param use_poly Logical. If TRUE, apply polynomial and interaction expansion to x_train before model fitting.
#'
#' @return A list containing the selected lambda value and the full cross-validated lasso model object.

tune_lasso <- function(x_train, y_train, family, use_poly = FALSE) {
  # transform features if needed
  if (use_poly) {
    x_train_poly <- transform_x_lasso(x_train)
  } else {
    x_train_poly <- x_train
  }
  
  # perform cross-validated lasso (alpha = 1)
  cv_lasso <- cv.glmnet(
    x = x_train_poly,
    y = y_train,
    alpha = 1,
    nfolds = 5,
    type.measure = "mse",
    standardize = TRUE,
    family = family
  )
  
  # select best lambda
  lambda <- cv_lasso$lambda.min
  
  # return selected lambda and model
  best_params <- list(
    lambda = lambda
  )
  
  return(best_params)
}
