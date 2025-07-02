BIC.drcOrdinal <- function(object, ...){
  dots <- list(...)
  if (!is.null(dots$epsilon)){
    epsilon <- dots$epsilon
  } else {
    epsilon <- 1e-16
  }
  
  n.parameters <- sum(sapply(object$drmList, function(mod) length(mod$coefficients)))
  n.obs <- sum(object$data[[object$weights]])
  n.parameters * log(n.obs) - 2 * logLik(object, epsilon) 
}
