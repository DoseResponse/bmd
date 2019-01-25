bmdBoot <- function(object, bmr, R=1000, boot="nonparametric", bmdType = "orig",
                    backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
                    backg=NA, 
                    def = c("excess", "additional", 
                            "relative", "extra", "added", "hybridExc", "hybridAdd", "point"),
                    bootInterval = c("percentile","BCa")){
  tmp.data <- bootDataGen(object,R,boot)
    drm.list <- lapply(tmp.data, function(x){
      drm(object$call$formula, data = x, type = object$type, fct = object[["fct"]])}
      )
      
    bmd.list <- lapply(drm.list,function(x){
    bmd(x, bmr = bmr, backgType = backgType, backg=backg, def=def, display=FALSE)[["Results"]][1]}
    )
    if(identical(object$type, "continuous")){
    if(identical(bootInterval,"BCa")){
      jackData <- list()
      for(i in 1:(dim(object$data)[1])){
        jackData[[i]] <- object$data[-i,]
      }
      bootJack <- sapply(jackData, function(x){
          bmd(drm(object$call$formula, data = x, type = object$type, fct = object[["fct"]]),
              bmr, backgType = backgType, backg = backg, def = def, interval = "delta", display=FALSE)$Results[1]
        }
        )
      use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, display=FALSE)[["Results"]][1]
      as.numeric(BCa(obs = use.bmd, data = object$data, unlist(bmd.list), bootJack)[1])
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
      bootJack <- sapply(jackData, function(x){
        bmd(drm(number~dose, data = x, type = "binomial", fct = object[["fct"]]),
            bmr, backgType = backgType, backg = backg, def = def, interval = "delta", display=FALSE)$Results[1]
      }
      )
      use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, display=FALSE)[["Results"]][1]
      BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = data.e, unlist(bmd.list), bootJack)[1])
    }
}
      if(bmdType == "orig"){
      use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, display=FALSE)[["Results"]][1]
    } else if(bmdType == "mean"){
      use.bmd <- mean(unlist(bmd.list))
    } else if(bmdType == "median"){
      use.bmd <- quantile(unlist(bmd.list),c(0.5))  
    } 
    
    BMDL <- ifelse(identical(bootInterval,"BCa"), BCaBMDL, quantile(unlist(bmd.list),c(0.05)))
    
    resMat <- matrix(NA,1,2)
    resMat[1,1] <- use.bmd
    resMat[1,2] <- BMDL
    colnames(resMat) <- c("BMD", "BMDL")
    rownames(resMat) <- c("")
    
    resBMD<-list(Results = resMat,
                 bootEst = unlist(bmd.list),
                 percentileInterval = quantile(unlist(bmd.list),c(0.05,0.95)))
    class(resBMD) = "bmd"
    return(resBMD)
   
}


