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