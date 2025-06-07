# check, install, and load necessary packages
packages <- c()
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))


#' Estimation of Treatment Effects and Standard Errors from Influence Functions Components
#'
#' Computes point estimates and standard errors for treatment effects based on arrays 
#' of influence function components from cross-fitting repetitions.
#' 
#' The function calculates the median and mean treatment effect across repetitions, 
#' as well as multiple standard error estimators (median-based, mean-based, 
#' with/without cross-fitting correction, influence function aggregation-based).
#'
#' @param psi_a A numeric array of shape (n_obs, n_rep, 1) containing the first influence function component.
#' @param psi_b A numeric array of shape (n_obs, n_rep, 1) containing the second influence function component.
#'
#' @return A named list containing:
#' \describe{
#'   \item{median_theta}{Median treatment effect estimate across repetitions.}
#'   \item{mean_theta}{Mean treatment effect estimate across repetitions.}
#'   \item{se_median}{Standard error based on median variance estimate (includes correction term).}
#'   \item{se_median_WO_CT}{Standard error based on median variance estimate (without correction term).}
#'   \item{se_mean}{Standard error based on mean variance estimate (includes correction term).}
#'   \item{se_mean_WO_CT}{Standard error based on mean variance estimate (without correction term).}
#'   \item{se_IF}{Standard error based on influence function aggregation.}
#' }
get_estimates <- function(psi_a, psi_b) {
  
  n_obs <- dim(psi_a)[1]
  n_rep <- dim(psi_a)[2]
  
  # sum over observations per repetition
  sum_psi_a <- apply(psi_a, 2, sum, na.rm = TRUE)
  sum_psi_b <- apply(psi_b, 2, sum, na.rm = TRUE)
  
  # estimate theta per repetition
  theta <- -(sum_psi_b) / (sum_psi_a)
  
  # point estimates
  median_theta <- median(theta, na.rm = TRUE)
  mean_theta <- mean(theta, na.rm = TRUE)
  
  # calculate Influence Function for each observation and repetition
  psi <- t(t(psi_a[,,1]) * theta + t(psi_b[,,1]))
  
  # variance of Influence Function for each repetition
  sigma <- apply(psi^2, 2, mean, na.rm = TRUE) / apply(psi_a, 2, mean, na.rm = TRUE)^2
  all_se <- sqrt(sigma / n_obs)  # standard error per repetition
  
  # variance across repetitions including cross-fitting variability
  var_se <- n_obs * all_se^2 + (theta - median(theta, na.rm = TRUE))^2
  var_se_WO_CT <- n_obs * all_se^2  
  
  # aggregated standard errors
  median_var <- median(var_se, na.rm = TRUE)
  se_median <- sqrt(median_var / n_obs)
  
  median_var_WO_CT <- median(var_se_WO_CT, na.rm = TRUE)
  se_median_WO_CT <- sqrt(median_var_WO_CT / n_obs)
  
  mean_var <- mean(var_se, na.rm = TRUE)
  se_mean <- sqrt(mean_var / n_obs)
  
  mean_var_WO_CT <- mean(var_se_WO_CT, na.rm = TRUE)
  se_mean_WO_CT <- sqrt(mean_var_WO_CT / n_obs)
  
  # Influence Function based standard errors
  Psi <- - sweep(psi, 2, apply(psi_a, 2, mean, na.rm = TRUE), "/")
  IF_combined <- apply(Psi, 1, mean, na.rm = TRUE)
  var <- var(IF_combined, na.rm = TRUE)
  se_IF <- sqrt(var / n_obs)
  
  return(list(
    median_theta = median_theta,
    mean_theta = mean_theta,
    se_median = se_median,
    se_median_WO_CT = se_median_WO_CT,
    se_mean = se_mean,
    se_mean_WO_CT = se_mean_WO_CT,
    se_IF = se_IF
  ))
}
