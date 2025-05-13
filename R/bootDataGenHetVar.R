bootDataGenHetVar <- function(object, R=1000, bootType=c("nonparametric", "parametric")){
  if(!inherits(object, "drcHetVar")){
    stop('bootData can only be generated from object of class "drcHetVar"')
  }
  bootType <- match.arg(bootType)
  
  tmp.data <- list()
  if(identical(bootType, "nonparametric")){
    dName <- object$dataList$names$dName
    data.e <- object$data
    data.e[,"row.num"] <- 1:nrow(data.e)
    data.e[,"dose"] <- data.e[,dName]
    tmp.data <- list()
    for(i in 1:R){
      tmp.data[[i]] <- data.e[as.numeric(unlist(aggregate(row.num ~ dose, data=data.e, 
                                                          FUN=function(x) sample(x,replace=TRUE))[[2]])),]
    }
  } else {
    dName <- object$dataList$names$dName
    rName <- object$dataList$names$rName
    
    data.e <- object$data
    data.e[,"dose"] <- data.e[,object$dataList$names$dName]
    for(i in 1:R){
      new.resp <- rnorm(nrow(data.e), mean = object$curve(data.e[,"dose"]), sd = object$sigmaFun(data.e[,"dose"]))
      data.e.copy <- data.e
      data.e.copy[,rName] <- new.resp
      tmp.data[[i]] <- data.e.copy
    }
  }
  
  return(tmp.data)
}