#' Log-Likelihood Method for drcHetVar Objects
#'
#' Extracts the log-likelihood from drcHetVar model objects.
#'
#' @param object A drcHetVar model object
#' @param ... Additional arguments (not used)
#'
#' @return The log-likelihood value
#' @export
logLik.drcHetVar <- function(object, ...){
  - object$sumList$numObs/2 * log(2*pi) - object$value/2
}

#' AIC Method for drcHetVar Objects
#'
#' Extracts the AIC from drcHetVar model objects.
#'
#' @param object A drcHetVar model object
#' @param ... Additional arguments (not used)
#' @param k Numeric value for penalty term in AIC calculation (default is 2)
#'
#' @return The log-likelihood value
#' @export
AIC.drcHetVar <- function(object, ..., k = 2){
  2*object$sumList$numPar - k*logLik(object)
}

#' BIC Method for drcHetVar Objects
#'
#' Extracts the BIC from drcHetVar model objects.
#'
#' @param object A drcHetVar model object
#' @param ... Additional arguments (not used)
#'
#' @return The log-likelihood value
#' @export
BIC.drcHetVar <- function(object, ...){
  object$sumList$numPar*log(object$sumList$numObs) - 2*logLik(object)
}
