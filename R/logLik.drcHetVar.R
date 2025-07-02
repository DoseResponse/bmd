logLik.drcHetVar <- function(object, ...){
  - object$sumList$numObs/2 * log(2*pi) - object$value/2
}

AIC.drcHetVar <- function(object, ..., k = 2){
  2*object$sumList$numPar - k*logLik(object)
}

BIC.drcHetVar <- function(object, ...){
  object$sumList$numPar*log(object$sumList$numObs) - 2*logLik(object)
}
