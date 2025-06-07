# check, install, and load necessary packages
packages <- c("dplyr", "data.table")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))


#' Variance Decomposition for Double ML Estimates
#'
#' Computes total variance and decomposes it into three components: variance across DGPs,
#' variance due to cross-fitting repetitions, and variance due to the learner itself.
#'
#' @param results_dml A data frame or data.table containing DML results with columns:
#'   - theta_hat: Estimated treatment effects
#'   - dgp: Identifier for the data generating process
#'   - cv_rep: Identifier for the cross-fitting repetition
#'   - rf_rep: Identifier for the learner repetition
#'
#' @return A data.table with:
#'   - variance_total: Total variance of theta_hat
#'   - variance_dgp: Between-DGP variance
#'   - share_variance_dgp: Share of total variance due to DGP differences
#'   - variance_cf_rep: Variance due to cross-fitting repetitions
#'   - share_variance_cf_rep: Share of total variance due to cross-fitting
#'   - variance_learner_rep: Residual learner-level variance
#'   - share_variance_learner_rep: Share of total variance due to learner noise
#'
variance_decomposition_sim1<- function(results_dml) {
  # total variance across all theta_hat estimates
  variance_total <- mean((results_dml$theta_hat - mean(results_dml$theta_hat))^2)
  
  # between-DGP variance: Variance of mean effects across different DGPs
  means_dgp <- results_dml %>%
    group_by(dgp) %>%
    summarise(mean_theta_hat_dgp = mean(theta_hat, na.rm = TRUE))
  variance_dgp <- mean((means_dgp$mean_theta_hat_dgp - mean(results_dml$theta_hat))^2)
  
  # variance due to cross-fitting repetitions (within each DGP)
  means_cross <- results_dml %>%
    group_by(dgp, cv_rep) %>%
    summarise(mean_theta_cross = mean(theta_hat, na.rm = TRUE))
  variances_cross <- means_cross %>%
    group_by(dgp) %>%
    summarise(variances_cross = mean((mean_theta_cross - mean(mean_theta_cross))^2))
  variance_cf_rep <- mean(variances_cross$variances_cross)
  
  # residual learner variance: Variation within each cross-fitting run
  variance_learner <- results_dml %>%
    group_by(dgp, cv_rep) %>%
    summarise(variance_learner = mean((theta_hat - mean(theta_hat))^2))
  variance_learner_rep <- mean(variance_learner$variance_learner)
  
  # normalize each component by the total variance
  share_variance_dgp <- variance_dgp / variance_total
  share_variance_cf_rep <- variance_cf_rep / variance_total
  share_variance_learner_rep <- variance_learner_rep / variance_total
  
  # return all variance components and shares
  return(data.table(
    variance_total = variance_total,
    variance_dgp = variance_dgp,
    share_variance_dgp = share_variance_dgp,
    variance_cf_rep = variance_cf_rep,
    share_variance_cf_rep = share_variance_cf_rep,
    variance_learner_rep = variance_learner_rep,
    share_variance_learner_rep = share_variance_learner_rep
  ))
}
