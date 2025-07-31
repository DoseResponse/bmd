bmd<-function(object, bmr, backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
               backg=NA, controlSD=NA,
               def = c("excess", "additional", 
                       "relative", "extra", "added", "hybridExc", "hybridAdd", "point"), 
               respTrans = c("none", "log", "sqrt"),
              interval = c("delta", "sandwich", "inv", "profile", "profileGrid"), sandwich.vcov=FALSE, display = TRUE, level=0.95, profileGridSize = NA, profileProgressInfo = TRUE) 
{
  if (missing(object)){
    stop(paste("object is missing", sep=""))
  } else {
    if(!inherits(object, "drc")){ stop('object must be of class "drc"')}
  }
  if (missing(def)) {
    stop(paste("def is missing", sep=""))
  }
  if(def=="point"){
    backgType <- "modelBased"
    } 
  if (missing(backgType)) {
    if(!(def %in% c("hybridExc", "hybridAdd"))){
      backgType <- "modelBased"
    } else {
      backgType <- "hybridSD"
    }
  }
  if (!(def %in% c("excess", "additional", "relative", "extra", "added", "hybridExc", "hybridAdd", "point"))) {
    stop(paste("Could not recognize def", sep=""))
  }
  if (!(backgType %in% c("modelBased","absolute","hybridSD","hybridPercentile"))) {
    stop(paste("Could not recognize backgType", sep=""))
  }
  
  level <- 1-2*(1-level)
  
  interval <- match.arg(interval)
  if(inherits(object, "drcMMRE") & interval != "delta"){
    stop("only delta type confidence interval supported for object of type \"drcMMRE\"")
  }
  if(interval == "sandwich"){
    sandwich.vcov <- TRUE
    interval <- "delta"
  }
  if(sandwich.vcov & !requireNamespace("sandwich")){
    stop('package "sandwich" must be installed to compute sandwich confidence intervals')
  }
  respTrans <- match.arg(respTrans)
  
  if(inherits(object$fct, "braincousens") & is.null(object$fct$fixed)){
    if(object$fct$name == "BC.4"){
      object$fct$fixed <- c(NA, 0, NA, NA, NA)
    } else if(object$fct$name == "BC.5"){
      object$fct$fixed <- c(NA, NA, NA, NA, NA)
    }
  }
  
  # Extract information from model
  # EDlist <- object$fct[["edfct"]] # Change after drc package has been updated with with edfct 
  EDlist <- bmd.edfct(object)
  parmMat <- object$parmMat
  nCurves <- ncol(parmMat)
  
  # bmrScaledList
  bmrScaledList <- getBmrScaledList(object, bmr, backgType, backg, controlSD, def, respTrans)
  
  fctDerivx <- object$fct$derivx
  if(is.null(fctDerivx)){
    fctDerivx <- getFctDerivx(object)
    if(is.null(fctDerivx) & (interval == "delta")) stop(paste0("Derivative of dose-response curve not defined for model: ", object$fct$name, "\nDelta confidence interval not available."))
  }
  
  # SINGLE CURVE
  if(nCurves == 1){
    
    # Initialise result matrices
    resMat <- matrix(NA, nrow = 1, ncol = 2, dimnames = list("", c("BMD", "BMDL")))
    bmrScaled <- matrix(bmrScaledList$bmrScaledMat, nrow = 1, ncol = 1, dimnames = list("", "bmrScaled"))
    bmdInterval <- matrix(NA, nrow = 1, ncol = 2, dimnames = list("", c("Lower", "Upper")))
    bmdSE <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("", "SE"))
    
    EDeval <- EDlist(parmMat[, 1], bmrScaled[1,], type = "absolute", reference = "control")
    bmdVal <- EDeval[[1]]
    bmdSEVal <- NA
    
    if(interval == "delta"){
      if(sandwich.vcov){
        varCov <- sandwich::sandwich(object)
      } else {
        varCov <- vcov(object)
      }
      dBmdVal <- EDeval[[2]] + bmrScaledList$dBmrScaled[,1] / fctDerivx(bmdVal, t(parmMat))[1]
      bmdSEVal <- sqrt(dBmdVal %*% varCov %*% dBmdVal)
      if(inherits(object, "drcMMRE")){
        intMat <- matrix(qnorm(c((1-level)/2, 1-(1-level)/2), mean = bmdVal, sd = bmdSEVal), ncol = 2)
      } else {
        intMat <- confint.basic(matrix(c(bmdVal, bmdSEVal), ncol = 2), 
                                      level = level, object$"type", df.residual(object), FALSE)
      }
    } else if(interval == "inv"){
      if(!identical(respTrans, "none")){stop("inverse regression interval not available for transformed response.")}
      slope <- drop(ifelse(object$curve[[1]](0)-object$curve[[1]](Inf)>0,"decreasing","increasing"))
      if(is.na(object$curve[[1]](0)-object$curve[[1]](Inf))){
        slope <- drop(ifelse(object$curve[[1]](0.00000001)-object$curve[[1]](100000000)>0,"decreasing","increasing"))
      }
      # useSD
      if(def %in% c("hybridAdd","hybridExc")){
        useSD <- ifelse(!is.na(controlSD),controlSD,sqrt(summary(object)$resVar))
      }
      
      intMat <- invBmd(object, bmr, level = level, slope=slope, backgType=backgType,
                       backg=backg, catLev=NA, extFactor=10, def=def, useSD=useSD,
                       sandwich.vcov=sandwich.vcov)[,c("BMDL", "BMDU"), drop = FALSE]
    } else if(interval == "profile"){
      tmpVals <- ED(object, bmrScaled, interval = "delta",
                        level = level, type = "absolute", vcov. = vcov, display = FALSE)[,c("Lower", "Upper"), drop = FALSE]
      if(backgType %in% c("modelBased", "absolute") & !(def %in% c("hybridExc", "hybridAdd"))
         & substr(object$fct$name,1,2) %in% c("LL", "LN", "W1", "W2")){
        if(is.na(profileGridSize)){
          profileGridSize <- 20
        }
        
        slope <- drop(ifelse(object$curve[[1]](0)-object$curve[[1]](Inf)>0,"decreasing","increasing"))
        if(is.na(object$curve[[1]](0)-object$curve[[1]](Inf))){
          slope <- drop(ifelse(object$curve[[1]](0.00000001)-object$curve[[1]](100000000)>0,"decreasing","increasing"))
        }
        
        tmpInterval <- bmdProfileCI(object, slope, bmr, backgType, backg, def, respTrans, level = level, gridSize = profileGridSize,
                                    bmdEst = bmdVal, lower = tmpVals[,"Lower"], upper = tmpVals[,"Upper"])
      } else {
        cat("\nReparametrised model not available for chosen backgType, def and dose-response curve. Proceeding with grid search.\n")
        if(missing(profileGridSize)){
          profileGridSize <- 50
        }
        tmpInterval <- bmdProfileCIgrid(object, bmr = bmr, backgType = backgType, def = def, level = level,
                                        gridSize = profileGridSize, progressInfo = profileProgressInfo)
      }
      intMat <- matrix(tmpInterval, ncol = 2)
    }
    
  if(interval == "profileGrid"){
    if(is.na(profileGridSize)){profileGridSize <- 50}
    tmpInterval <- bmdProfileCIgrid(object, bmr = bmr, backgType = backgType, def = def, level = level,
                                    gridSize = profileGridSize, progressInfo = profileProgressInfo)
    intMat <- matrix(tmpInterval, ncol = 2)
  }
    
    resMat[1,] <- c(bmdVal, intMat[1,1])
    bmdInterval[1,] <- intMat
    bmdSE[1,] <- bmdSEVal
  }
  
  # MULTIPLE CURVES FITTED
  else {
    if(is.null(object$objList)){
      # CURVES ARE NOT FITTED INDEPENDENTLY
      curveNames <- colnames(parmMat)
      if(sandwich.vcov){
        vcMat <- sandwich::sandwich(object)
      } else {
        vcMat <- vcov(object)
      }
      
      # Initialise result matrices
      resMat <- matrix(NA, nrow = nCurves, ncol = 2, dimnames = list(curveNames, c("BMD", "BMDL")))
      bmrScaled <- matrix(bmrScaledList$bmrScaledMat, nrow = nCurves, ncol = 1, dimnames = list(curveNames, "bmrScaled"))
      bmdInterval <- matrix(NA, nrow = nCurves, ncol = 2, dimnames = list(curveNames, c("Lower", "Upper")))
      bmdSE <- matrix(NA, nrow = nCurves, ncol = 1, dimnames = list(curveNames, "SE"))
      
      # pmodelsMatrixList - used for constructing correct vcov matrices
      pmodelsMatrixList <- getpmodelsMatrixList(object)
      
      for(iCurve in 1:nCurves){
        parmChosen <- parmMat[, iCurve]
        varCov <- pmodelsMatrixList[[iCurve]] %*% vcMat %*% t(pmodelsMatrixList[[iCurve]])
        
        EDeval <- EDlist(parmChosen, bmrScaled[iCurve,], type = "absolute", reference = "control") 
        bmdVal <- EDeval[[1]]
        dBmdVal <- EDeval[[2]] + bmrScaledList$dBmrScaled[,iCurve] / fctDerivx(bmdVal, t(parmMat))[iCurve]
        bmdSEVal <- sqrt(dBmdVal %*% varCov %*% dBmdVal)
        
        if(interval == "delta"){
          intMat <- confint.basic(matrix(c(bmdVal, bmdSEVal), ncol = 2), 
                                        level = level, object$"type", df.residual(object), FALSE)
        }
        
        resMat[iCurve,] <- c(bmdVal, intMat[1,1])
        bmdInterval[iCurve,] <- intMat
        bmdSE[iCurve,] <- bmdSEVal
      }
    } else {
      # CURVES ARE FITTED INDEPENDENTLY
      bmdCall <- function(object){
        bmd(object, bmr, backgType, backg, controlSD,
            def, respTrans, interval, sandwich.vcov, display = FALSE, level = 1 + (level-1)/2, profileGridSize, profileProgressInfo)
      }
      
      bmdList <- lapply(object$objList, bmdCall)
      
      resMat <- do.call(rbind, lapply(bmdList, function(x) x$Results))
      bmrScaled <- do.call(rbind, lapply(bmdList, function(x) x$bmrScaled))
      bmdInterval <- do.call(rbind, lapply(bmdList, function(x) x$interval))
      bmdSE <- do.call(rbind, lapply(bmdList, function(x) x$SE))
      
      # set rownames of result matrices
      rownames(resMat) <- levels(object$dataList$curveid)
      rownames(bmrScaled) <- levels(object$dataList$curveid)
      rownames(bmdInterval) <- levels(object$dataList$curveid)
      rownames(bmdSE) <- levels(object$dataList$curveid)
    }
  }
  
  if (display) {
    print(resMat)
  }
  
  resBMD<-list(Results = resMat,
               bmrScaled = bmrScaled,
               interval = bmdInterval,
               SE = bmdSE,
               model = object)
  class(resBMD) <- "bmd"
  invisible(resBMD) 
}



