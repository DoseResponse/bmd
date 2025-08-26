
#' AIC Method for drcOrdinal Objects
#'
#' Calculates the Akaike Information Criterion for drcOrdinal model objects.
#'
#' @param object A drcOrdinal model object
#' @param ... Additional arguments (not used)
#' @param k The penalty per parameter to be used; default is 2
#'
#' @return The AIC value
#' @export
AIC.drcOrdinal <- function(object, ..., k = 2) {
  dots <- list(...)
  if (!is.null(dots$epsilon)){
    epsilon <- dots$epsilon
  } else {
    epsilon <- 1e-16
  }
  
  n.parameters <- sum(sapply(object$drmList, function(mod) length(mod$coefficients)))
  -2 * logLik(object, epsilon = epsilon) + k * n.parameters
}
