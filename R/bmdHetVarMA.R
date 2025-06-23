bmdHetVarMA <- function(modelList, modelWeights = c("AIC", "BIC"), bmr, backgType = c("absolute", "hybridSD", "hybridPercentile"), backg = NA, def = c("hybridExc", "hybridAdd"), interval = c("boot", "none"), R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE){
  ### Assertions ###
  # modelList
  if(any(!sapply(modelList, inherits,"drcHetVar"))){
    stop('modelList must be a list of dose-response models with a heterogeneous variance structure of class "drcHetVar" ')
  }
  
  # modelWeights
  if(!(length(modelWeights) %in% c(1, length(modelList)))){
    stop('modelWeights must either be "AIC", "BIC" or a numeric vector of same length as modelList')
  } else if(length(modelWeights) == 1){
    if(!(modelWeights %in% c("AIC", "BIC"))) stop('modelWeights must either be "AIC", "BIC" or a numeric vector of same length as modelList')
  } else if(any(!is.numeric(modelWeights))){
    stop('modelWeights must either be "AIC", "BIC" or a numeric vector of same length as modelList')
  }
  
  # bmr
  if(missing(bmr)){
    stop('argument "bmr" needs to be specified as a number between 0 and 1')
  }
  if(!is.numeric(bmr)){
    stop('argument "bmr" needs to be specified as a number between 0 and 1')
  }
  if(bmr <= 0 | bmr >=1){
    stop('argument "bmr" needs to be specified as a number between 0 and 1')
  }
  
  # backgType
  if (missing(backgType)) {
    stop('backgType is missing. Options are "absolute", "hybridSD" or "hybridPercentile"')
  }
  if (!(backgType %in% c("absolute","hybridSD","hybridPercentile"))) {
    stop('Could not recognize backgType. Options are "absolute", "hybridSD" or "hybridPercentile"')
  }
  
  # def
  if(missing(def)){
    stop('def is missing. Options are "hybridExc" or "hybridAdd"')
  }
  if(!def %in% c("hybridExc", "hybridAdd")){
    stop('Could not recognize def. Options are "hybridExc" or "hybridAdd"')
  }
  
  level <- 1-2*(1-level)
  
  #bmdHetVarList
  bmdHetVarList <- lapply(modelList, 
                          FUN=function(object){
                            bmdHetVar(object, bmr = bmr, backgType = backgType, backg = backg, def = def,
                                                      interval = "none", display=FALSE)})
  if(any(sapply(bmdHetVarList, inherits, "try-error"))){
    fail_models <- which(sapply(bmdHetVarList, inherits, "try-error"))
    stop(paste0("bmd could not be estimated for the following model(s): ", paste0(fail_models, collapse = ", ")))
  }
  
  # Estimate weights
  if(identical(modelWeights,"AIC")){
    modelWeights0<-exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))/2)/
      sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))/2))
  } else if(identical(modelWeights,"BIC")){
    modelWeights0<-exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))/2)/
      sum(exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))/2))
  } else {
    modelWeights0 <- modelWeights/sum(modelWeights)
  }
  
  # bmdEst
  bmdVals <- sapply(bmdHetVarList, function(z) z$Results[1])
  bmdEst <- sum(bmdVals * modelWeights0)
  
  # INTERVAL
  interval <- match.arg(interval)
  if(identical(interval, "none")){
    BMDL <- NA
    BMDU <- NA
    used.Boot <- NA
  } else {
    bootDataList <- bootDataGenHetVar(modelList[[1]], R = R, bootType = "nonparametric")
    
    bmdHetVarMABoot <- function(bootData){
      bootModHetVarList <- lapply(modelList, function(object) drmHetVar(formula = object$formula, var.formula = object$var.formula, data = bootData, fct = object$fct))
      bootBmdEst <- bmdHetVarMA(modelList = bootModHetVarList, modelWeights = modelWeights, bmr = bmr, backgType = backgType, backg = backg,
                              def = def, interval = "none", display = FALSE)$Results[,1]
      bootBmdEst
    }
    
    if(progressInfo){
      cat("Performing bootstrap\n")
      pb <- txtProgressBar(min = 0, max = R, style = 3)
    }
    
    bootBmdEst <- numeric(R)
    for(i in 1:R){
      bootBmdEst[i] <- suppressWarnings(as.numeric(try(bmdHetVarMABoot(bootDataList[[i]]), silent = TRUE)))
      if(progressInfo) setTxtProgressBar(pb, i)
    }
    if(progressInfo) close(pb)
    
    boot0<-bootBmdEst[!is.na(bootBmdEst)]
    used.Boot <- length(boot0)
    
    if(length(boot0) == 0){ 
      BMDL <- NA 
      BMDU <- NA
    } else {
      BMDL <- quantile(boot0,p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
      BMDU <- quantile(boot0,p=c(level), na.rm = TRUE)
    }
  }
  
  resMat <- matrix(c(bmdEst, BMDL), nrow = 1, ncol = 2, dimnames = list(NULL, c("BMD_MA", "BMDL_MA")))
  bmdInterval <- matrix(c(BMDL, BMDU), nrow = 1, ncol = 2, dimnames = list("", c("Lower", "Upper")))
  
  # GATHER RESULTS
  if (display) {
    print(resMat)
  }
  
  resBMD<-list(Results = resMat,
               Boot.samples.used = used.Boot,
               interval = bmdInterval,
               modelWeights = modelWeights0)
  class(resBMD) <- c("bmdHetVar", "bmd")
  invisible(resBMD) 
}