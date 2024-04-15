bmdBoot <- function(object, bmr, R=1000, bootType="nonparametric", bmdType = "orig",
                    backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
                    backg=NA, 
                    controlSD=NA,
                    def = c("excess", "additional", 
                            "relative", "extra", "added", "hybridExc", "hybridAdd", "point"),
                    respTrans = c("none", "log", "sqrt"),
                    bootInterval = c("percentile","BCa"),
                    display=TRUE, level=0.95){
  if (identical(object$type,"binomial") & bootType=="semiparametric") {
    stop(paste("\"Semiparametric bootstrap does not work for quantal data\"", sep=""))
  }
  if (object$type %in% c("Poisson","negbin1","negbin2") & bootType!="nonparametric") {
    stop(paste("\"",object$type,"\" only works with nonparametric bootstrap\"", sep=""))
  }
 
  if (object$type %in% c("binomial","continuous")) {
  
  tmp.data <- bootDataGen(object,R,bootType,aggregated=FALSE)
  
  
  drm.list.tmp <- lapply(tmp.data, function(x){
    try(drm(object$call$formula, data = x, type = object$type, fct = object[["fct"]]),TRUE)}
  )
  list.condition <- sapply(drm.list.tmp, function(x) class(x)=="drc")
  drm.list  <- drm.list.tmp[list.condition]
  
  bmd.list.try <- lapply(drm.list,function(x){
    try(bmd(x, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
            display=FALSE, level=level)[["Results"]][1],TRUE)}
  )
  bmd.list <- bmd.list.try[sapply(bmd.list.try, function(x) class(x)=="numeric")]
  }
  
  if (object$type %in% c("Poisson","negbin1","negbin2")) {
    tmp.data <- bootDataGen(object,R,bootType,aggregated=FALSE)
  
    drm.list.tmp <- lapply(tmp.data, function(x){
      try(drm(object$call$formula, data = x, type = object$type, weights=weights, fct = object[["fct"]]),TRUE)}
    )
    list.condition <- sapply(drm.list.tmp, function(x) class(x)=="drc")
    drm.list  <- drm.list.tmp[list.condition]
    
    bmd.list <- lapply(drm.list,function(x){
      bmd(x, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
          display=FALSE, level=level)[["Results"]][1]}
    )
    }
  
  
  
  if(identical(object$type, "continuous")){
    if(identical(bootInterval,"BCa")){
      jackData <- list()
      for(i in 1:(dim(object$data)[1])){
        jackData[[i]] <- object$data[-i,]
      }
      bootJack.drm.tmp <- lapply(jackData, function(x){
        try(drm(object$call$formula, data = x, fct = object[["fct"]]),TRUE)
      })
      list.condition <- sapply(bootJack.drm.tmp, function(x) class(x)=="drc")
      bootJack.drm<- bootJack.drm.tmp[list.condition]
      
      bootJack <- sapply(bootJack.drm, function(x){
        bmd(x, bmr, backgType = backgType, backg = backg, def = def, controlSD=controlSD, respTrans = respTrans,
            interval = "delta", display=FALSE, level=level)$Results[1]
      }
      )
      
      use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
                     display=FALSE, level=level)[["Results"]][1]
      BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = object$data, unlist(bmd.list), bootJack, level=level)[1])
      
    }
  }
  if(identical(object$type, "binomial")){
    if(identical(bootInterval, "BCa")){
      data.str <- object$data
      data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
      data.e<-expandBinomial(data.str, 
                             number = "number",
                             total = "weights",
                             dose = as.character(object$call$formula[[3]]))
      jackData <- list()
      for(i in 1:(dim(data.e)[1])){
        jackData[[i]] <- data.e[-i,]
      }
      bootJack.drm.tmp <- lapply(jackData, function(x){
        try(drm(number~dose, data = x, type = "binomial", fct = object[["fct"]]),TRUE)
      })
      list.condition <- sapply(bootJack.drm.tmp, function(x) class(x)=="drc")
      bootJack.drm<- bootJack.drm.tmp[list.condition]
      
      bootJack <- sapply(bootJack.drm, function(x){
        bmd(x, bmr, backgType = backgType, backg = backg, def = def, controlSD=controlSD, respTrans = respTrans,
            interval = "delta", display=FALSE, level=level)$Results[1]
      }
      )
      
      use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
                     display=FALSE, level=level)[["Results"]][1]
      BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = data.e, unlist(bmd.list), bootJack, level=level)[1])
    }
  }
  if(object$type %in% c("Poisson","negbin1","negbin2")){
    if(identical(bootInterval,"BCa")){
      jackData <- list()
      for(i in 1:(dim(object$data)[1])){
        jackData[[i]] <- object$data[-i,]
      }
      bootJack.drm.tmp <- lapply(jackData, function(x){
        try(drm(object$call$formula, data = x, type = object$type, weights=weights, fct = object[["fct"]]),TRUE)
      })
      list.condition <- sapply(bootJack.drm.tmp, function(x) class(x)=="drc")
      bootJack.drm<- bootJack.drm.tmp[list.condition]
      
      bootJack <- sapply(bootJack.drm, function(x){
        bmd(x, bmr, backgType = backgType, backg = backg, def = def, controlSD=controlSD, interval = "delta", respTrans = respTrans,
            display=FALSE, level=level)$Results[1]
      }
      )
      
      use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
                     display=FALSE, level=level)[["Results"]][1]
      BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = object$data, bootSample=unlist(bmd.list), 
                                bootjack=bootJack, level=level)[1])
    }
  }
  if(bmdType == "orig"){
    use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
                   display=FALSE, level=level)[["Results"]][1]
  } else if(bmdType == "mean"){
    use.bmd <- mean(unlist(bmd.list))
  } else if(bmdType == "median"){
    use.bmd <- quantile(unlist(bmd.list),c(1-level))  
  } 
  
  BMDL <- ifelse(identical(bootInterval,"BCa"), BCaBMDL, quantile(unlist(bmd.list),c(1-level),na.rm = TRUE))
  BMDU <- ifelse(identical(bootInterval,"BCa"), "Not available for BCa bootstrap", quantile(unlist(bmd.list),c(level),na.rm = TRUE))
  
  resMat <- matrix(NA,1,2)
  resMat[1,1] <- use.bmd
  resMat[1,2] <- BMDL
  colnames(resMat) <- c("BMD", "BMDL")
  rownames(resMat) <- c("")
  
  if(display){
    print(resMat)
  }
  
  resBMD<-list(Results = resMat,
               bootEst = unlist(bmd.list),
               Interval = c(BMDL, BMDU))
  class(resBMD) <- "bmd"
  invisible(resBMD)
}