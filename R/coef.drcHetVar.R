coef.drcHetVar <- function(object){
  curvePar <- object$curvePar
  names(curvePar) <- paste0("Curve:", names(curvePar))
  sigmaPar <- object$sigmaPar
  names(sigmaPar) <- paste0("Sigma:", names(sigmaPar))
  c(curvePar, sigmaPar)
}