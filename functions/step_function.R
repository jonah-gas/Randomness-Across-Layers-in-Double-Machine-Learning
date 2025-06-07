# Check, install, and load necessary packages
packages <- c()
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(sapply(packages, install_if_missing))
invisible(lapply(packages, library, character.only = TRUE))

#' Step function to induce non-linearity in a numeric vector
#'
#' @param x numeric vector. Input values to be transformed.
#'
#' @return A numeric vector of the same length with transformed values:
#'   - -1 if value ≤ 25% quantile
#'   -  0 if value ≤ 50% quantile
#'   -  0.5 if value ≤ 75% quantile
#'   -  1 otherwise
#'
step_function <- function(x) {
  # calculate quartiles of the input vector
  quartiles <- quantile(x, probs = c(0.25, 0.5, 0.75))
  
  # return step values based on quartile cutoffs
  sapply(x, function(val) {
    if (val <= quartiles[1]) {
      return(-1)
    } else if (val <= quartiles[2]) {
      return(0)
    } else if (val <= quartiles[3]) {
      return(0.5)
    } else {
      return(1)
    }
  })
}
