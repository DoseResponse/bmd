logLik.drcOrdinal <- function(object, epsilon = 10^(-16)){
  tmp <- sapply(1:length(object$levels), function(cat.i){
    cat <- object$levels[[cat.i]]
    cat.per.dose <- unlist(object$drmList[[1]]$origData[,cat])
    if(cat.i == 1){
      p.fun <- function(x) 1 - object$drmList[[cat.i]]$curve[[1]](x)
    } else if(cat.i == length(object$levels)){
      p.fun <- function(x) object$drmList[[cat.i-1]]$curve[[1]](x)
    } else{
      p.fun <- function(x) object$drmList[[cat.i-1]]$curve[[1]](x) - object$drmList[[cat.i]]$curve[[1]](x)
    }
    p.val <- sapply(object$drmList[[1]]$origData[, as.character(object$drmList[[1]]$call$formula[[3]])],
                    function(x) pmax(epsilon, p.fun(x)))
    log(p.val) * cat.per.dose
  }
  )
  sum(tmp)
}
