# check, install, and load necessary packages
packages <- c()
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))


#' DML Inference Function for Simulation 1
#'
#' Computes point estimate, standard error, t-statistic, and p-value based on orthogonal score components.
#'
#' @param psi_a numeric vector. Score component (derivative) with respect to the target parameter theta.
#' @param psi_b numeric vector. Score component (residual) orthogonal to theta.
#'
#' @return A named numeric vector containing:
#'   - theta: Point estimate of the treatment effect
#'   - se: Standard error of the estimate
#'   - t: t-statistic for the hypothesis test
#'   - p: p-value associated with the test statistic
#'
dml_inference_sim1 <- function(psi_a, psi_b) {
  N <- length(psi_a)
  theta <- -sum(psi_b) / sum(psi_a)
  psi <- theta * psi_a + psi_b
  Psi <- -psi / mean(psi_a)
  sigma2 <- var(Psi)
  se <- sqrt(sigma2 / N)
  t <- theta / se
  p <- 2 * pt(abs(t), df = N, lower.tail = FALSE)
  
  return(c(
    theta = theta,
    se = se,
    t = t,
    p = p
  ))
}
