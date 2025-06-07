# check, install, and load necessary packages
packages <- c()
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))

#' Data Generating Process 1 (DGP1) for treatment effect simulation
#'
#' @param n_obs integer. Number of observations to generate.
#' @param theta numeric. True average treatment effect (ATE) used in the DGP.
#'
#' @return A named list containing:
#'   - y: Observed outcomes (numeric vector)
#'   - w: Treatment assignment (binary vector)
#'   - x: Covariates (matrix of size n_obs × 5)
#'   - theta: True treatment effect
#'   - n_obs: Number of observations
#'   - w_prob: Treatment probabilities (propensity scores)
#'
generate_DGP_1 <- function(n_obs, theta = 0.5 ) {
  n_p <- 5
  se_of_errorterm <- 1
  
  # generate covariates: U(-π, π)
  x <- matrix(runif(n_obs * n_p, -pi, pi), ncol = n_p)
  
  # define nuisance functions and treatment effect
  e <- function(x) pnorm(sin(x))              # propensity score
  m0 <- function(x) sin(x)                    # outcome under control
  m1 <- function(x) m0(x) + theta             # outcome under treatment
  
  # treatment assignment
  w_prob <- e(x[, 1])
  w <- rbinom(n_obs, 1, w_prob)
  
  # observed outcome
  y <- w * m1(x[, 1]) + (1 - w) * m0(x[, 1]) + rnorm(n_obs, 0, se_of_errorterm)
  
  # return simulated data
  return(list(
    y = y,
    w = w,
    x = x,
    theta = theta,
    n_obs = n_obs,
    w_prob = w_prob
  ))
}
