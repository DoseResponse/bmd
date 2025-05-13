bmdHetVar <- function(object, bmr, backgType = c("absolute", "hybridSD", "hybridPercentile"), backg = NA, def = c("hybridExc", "hybridAdd"), interval = c("boot", "none"), R = 1000, level = 0.95, bootType = "nonparametric", progressInfo = TRUE, display = TRUE){
  ### Assertions ###
  # object
  if(!inherits(object,"drcHetVar")){
    stop('object must be a dose-response model with a heterogeneous variance structure of class "drcHetVar" ')
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
  
  # SLOPE
  slope <- drop(ifelse(object$curve(0)-object$curve(Inf)>0,"decreasing","increasing"))
  if(is.na(object$curve(0)-object$curve(Inf))){
    slope <- drop(ifelse(object$curve(0.00000001)-object$curve(100000000)>0,"decreasing","increasing"))
  }
  
  # sigmaFun
  sigmaFun0 <- object$sigmaFun # sigmaFun(object, var.formula)
  
  # bmrScaled
  if(slope == "increasing"){
    # BACKGROUND
    if (identical(backgType,"absolute")) {
      if(is.na(backg)){
        stop('backgType = absolute, but backg not supplied')
      }
      p0 <- 1 - pnorm((backg - object$curve(0)) / sigmaFun0(0))
    }
    if(identical(backgType, "hybridPercentile")) {
      p0 <- ifelse(is.na(backg),1-0.9,1-backg)
    }
    if (identical(backgType,"hybridSD")) {
      p0 <- ifelse(is.na(backg), 1-pnorm(2), 1-pnorm(backg))
    }
    
    # BMRSCALED
    bmrScaled <- switch(
      def,
      hybridExc = function(x){ sigmaFun0(x) * 
          (qnorm(1 - p0) - qnorm(1 - p0 - (1 - p0)*bmr)) + object$curve(0)},
      hybridAdd = function(x){ sigmaFun0(x) * 
          (qnorm(1 - p0) - qnorm(1 - (p0 + bmr))) + object$curve(0)}
    ) 
  } else {
    # BACKGROUND
    if (identical(backgType,"absolute")) {
      if(is.na(backg)){
        stop('backgType = absolute, but backg not supplied')
      }
      p0 <- pnorm((backg - object$curve(0)) / sigmaFun0(0))
    }
    if(identical(backgType, "hybridPercentile")) {
      p0 <- ifelse(is.na(backg),0.1,backg)
    }
    if (identical(backgType,"hybridSD")) {
      p0 <- ifelse(is.na(backg), pnorm(-2), pnorm(-backg))
    }
    
    # BMRSCALED
    bmrScaled <- switch(
      def,
      hybridExc = function(x){ sigmaFun0(x) * 
          (qnorm(p0) - qnorm(bmr + (1-bmr) * p0)) + object$curve(0)},
      hybridAdd = function(x){ sigmaFun0(x) * 
          (qnorm(p0) - qnorm(bmr + p0)) + object$curve(0)}
    ) 
  }
  
  # BMD ESTIMATION
  f0 <- function(x) object$curve(x) - bmrScaled(x)
  interval0 <- range(object$dataList$dose, na.rm = TRUE)
  uniroot0 <- try(uniroot(f = f0, interval = interval0), silent = TRUE)
  
  if(inherits(uniroot0, "try-error")){
    bmdEst <- NA
    warning('error when estimating bmd. Root not found.\n')
  } else {
    bmdEst <- uniroot0$root
  }
  
  # INTERVAL
  interval <- match.arg(interval)
  if(identical(interval, "none")){
    BMDL <- NA
    BMDU <- NA
  } else {
    # drc_obj <- eval(substitute(drm(formula0, data = object$data, fct = fct0, type = "continuous", control = drmc(maxIt = 1, noMessage = TRUE)),
    #                            list(formula0 = object$formula,
    #                                 fct0 = object$fct
    #                            )))
    # bootDataList <- bootDataGen(drc_obj, R=R, bootType="nonparametric",aggregated=FALSE)
    bootDataList <- bootDataGenHetVar(object, R = R, bootType = bootType)
    
    bmdHetVarBoot <- function(bootData){
      bootModHetVar <- drmHetVar(formula = object$formula, var.formula = object$var.formula, data = bootData, fct = object$fct)
      bootBmdEst <- bmdHetVar(object = bootModHetVar, bmr = bmr, backgType = backgType, backg = backg,
                              def = def, interval = "none", display = FALSE)$Results[,1]
      bootBmdEst
    }
    
    if(progressInfo){
      cat("Performing bootstrap\n")
      pb <- txtProgressBar(min = 0, max = R, style = 3)
    }
    
    bootBmdEst <- numeric(R)
    for(i in 1:R){
      bootBmdEst[i] <- suppressWarnings(as.numeric(try(bmdHetVarBoot(bootDataList[[i]]), silent = TRUE)))
      if(progressInfo) setTxtProgressBar(pb, i)
    }
    if(progressInfo) close(pb)
    
    boot0<-bootBmdEst[!is.na(bootBmdEst)]
    
    if(length(boot0) == 0){ 
      BMDL <- NA 
      BMDU <- NA
    } else {
      BMDL <- quantile(boot0,p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
      BMDU <- quantile(boot0,p=c(level), na.rm = TRUE)
    }
  }
  
  resMat <- matrix(c(bmdEst, BMDL), nrow = 1, ncol = 2, dimnames = list(NULL, c("BMD", "BMDL")))
  bmrScaled <- matrix(object$curve(bmdEst), nrow = 1, ncol = 1, dimnames = list("", "bmrScaled"))
  bmdInterval <- matrix(c(BMDL, BMDU), nrow = 1, ncol = 2, dimnames = list("", c("Lower", "Upper")))
  
  # GATHER RESULTS
  if (display) {
    print(resMat)
  }
  
  resBMD<-list(Results = resMat,
               bmrScaled = bmrScaled,
               interval = bmdInterval,
               model = object)
  class(resBMD) <- "bmdHetVar"
  invisible(resBMD) 
}