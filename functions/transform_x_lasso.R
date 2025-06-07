# check, install, and load necessary packages
packages <- c()
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))

#' Polynomial and Pairwise Interaction Feature Expansion for Lasso
#'
#' Generates polynomial features and pairwise interaction terms for a given feature matrix.
#' Designed for use with lasso regression to capture nonlinearities and two-way interactions.
#'
#' @param x A numeric matrix or data.frame of input features.
#' @param powers A numeric vector indicating the polynomial degrees to include (default is 1:7).
#'
#' @return A matrix containing the expanded feature set, including polynomial terms and pairwise interactions.
transform_x_lasso <- function(x, powers = 1:7) {
  # polynomials
  x_poly <- sapply(powers, function(p) x^p, simplify = "array")
  x_poly <- do.call(cbind, lapply(seq_along(powers), function(i) x_poly[,,i]))
  
  # interactions
  interaction_terms <- apply(combn(ncol(x), 2), 2, function(idx) x[, idx[1]] * x[, idx[2]])
  
  # combine polynomial and interactions
  x_expanded <- cbind(x_poly, interaction_terms)
  
  # return expanded feature set
  return(x_expanded)
}
