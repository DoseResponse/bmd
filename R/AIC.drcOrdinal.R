AIC.drcOrdinal <- function(object, epsilon = 10^(-16)){
  n.parameters <- sum(sapply(object$drmList, function(mod) length(mod$coefficients)))
  - 2 * logLik(object, epsilon) + 2 * n.parameters
}
