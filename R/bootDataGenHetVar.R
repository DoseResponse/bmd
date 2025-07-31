bootDataGenHetVar <- function(object, R=1000, bootType=c("nonparametric", "semiparametric", "parametric")){
  if(!inherits(object, "drcHetVar")){
    stop('bootData can only be generated from object of class "drcHetVar"')
  }
  bootType <- match.arg(bootType)
  
  tmp.data <- list()
  if(identical(bootType, "nonparametric")){
    dName <- object$dataList$names$dName
    data.e <- as.data.frame(object$data)
    row.num <- 1:nrow(data.e)
    dose <- data.e[,dName]
    tmp.data <- list()
    for(i in 1:R){
      tmp.data[[i]] <- data.e[as.numeric(unlist(aggregate(row.num ~ dose, # data=data.e, 
                                                          FUN=function(x) sample(x,replace=TRUE))[[2]])),]
    }
  } else if(identical(bootType, "semiparametric")){
    dName <- object$dataList$names$dName
    rName <- object$dataList$names$rName
    
    data.e <- as.data.frame(object$data)
    dose <- data.e[,dName]
    sigma_fitted <- object$sigmaFun(dose)
    std_residuals <- object$residuals/sigma_fitted
    for(i in 1:R){
      new.std_residuals <- sample(std_residuals, replace = TRUE)
      new.resp <- object$curve(dose) + new.std_residuals * sigma_fitted
      data.e.copy <- data.e
      data.e.copy[,rName] <- new.resp
      tmp.data[[i]] <- data.e.copy
    }
  } else if(identical(bootType, "parametric")){
    dName <- object$dataList$names$dName
    rName <- object$dataList$names$rName
    
    data.e <- as.data.frame(object$data)
    dose <- data.e[,object$dataList$names$dName]
    fitted <- object$fitted
    sigma_fitted <- object$sigmaFun(dose)
    for(i in 1:R){
      new.resp <- rnorm(nrow(data.e), mean = fitted, sd = sigma_fitted)
      data.e.copy <- data.e
      data.e.copy[,rName] <- new.resp
      tmp.data[[i]] <- data.e.copy
    }
  }
  
  return(tmp.data)
}
