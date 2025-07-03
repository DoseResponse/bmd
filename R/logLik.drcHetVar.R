
#' @title S3 method
#' @export
logLik.drcHetVar <- function(object, ...){
  - object$sumList$numObs/2 * log(2*pi) - object$value/2
}

#' @title S3 method
#' @export
AIC.drcHetVar <- function(object, ..., k = 2){
  2*object$sumList$numPar - k*logLik(object)
}

#' @title S3 method
#' @export
BIC.drcHetVar <- function(object, ...){
  object$sumList$numPar*log(object$sumList$numObs) - 2*logLik(object)
}
