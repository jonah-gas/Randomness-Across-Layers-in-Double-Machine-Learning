# check, install, and load necessary packages
packages <- c("MASS")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))

#' Data Generating Process 2 (DGP2) for treatment effect simulation
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
generate_DGP_2 <- function(n_obs, theta = 1) {
  # number of covariates
  n_p <- 5
  
  # coefficient vectors for nuisance functions
  coef_e <- rep(0.1, n_p)
  coef_m0 <- rep(0.1, n_p)
  
  # random intercepts for propensity and outcome models
  alpha_0 <- rnorm(1, mean = 0, sd = 1)
  alpha_1 <- rnorm(1, mean = 0, sd = 1)
  
  # error terms for treated and untreated groups
  epsilon_0 <- rnorm(n_obs, mean = 0, sd = 1)
  epsilon_1 <- rnorm(n_obs, mean = 0, sd = 1)
  
  # generate covariates with multivariate normal distribution
  A <- matrix(rnorm(n_p * n_p, mean = 0, sd = sqrt(0.5)), nrow = n_p)
  Sigma <- t(A) %*% A
  x <- MASS::mvrnorm(n_obs, mu = rep(0, n_p), Sigma = Sigma)
  
  # apply step function to x3 to introduce non-linearity
  x[, 3] <- step_function(x[, 3])
  
  # define nuisance functions
  e <- function(x, coef_e, alpha_0) {
    pnorm(
      alpha_0 +
        coef_e[1] * x[, 1] +
        coef_e[2] * (x[, 2]^2) +
        coef_e[3] * (x[, 1] * x[, 2]) +
        coef_e[4] * x[, 3] +
        coef_e[5] * (x[, 4]^3)
    )
  }
  
  m0 <- function(x, coef_m0, alpha_1) {
    alpha_1 +
      coef_m0[1] * x[, 1] +
      coef_m0[2] * (x[, 2]^2) +
      coef_m0[3] * (x[, 1] * x[, 2]) +
      coef_m0[4] * x[, 3] +
      coef_m0[5] * (x[, 4]^3)
  }
  
  m1 <- function(x, coef_m0, alpha_1, theta) {
    m0(x, coef_m0, alpha_1) + theta
  }
  
  # treatment assignment
  w_prob <- e(x, coef_e, alpha_0)
  w_prob <- pmax(pmin(w_prob, 0.9), 0.1)  # truncation for stability
  w <- rbinom(n_obs, size = 1, prob = w_prob)
  
  # observed outcome
  y <- w * m1(x, coef_m0, alpha_1, theta) +
    (1 - w) * m0(x, coef_m0, alpha_1) +
    w * epsilon_1 + (1 - w) * epsilon_0
  
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
