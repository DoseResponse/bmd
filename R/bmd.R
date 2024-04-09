bmd<-function (object, bmr, backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
               backg=NA, controlSD=NA,
               def = c("excess", "additional", 
                       "relative", "extra", "added", "hybridExc", "hybridAdd", "point"), 
              interval = c("delta", "inv", "profile", "profileGrid"), sandwich.vcov=FALSE, display = TRUE, level=0.95, profileGridSize, profileProgressInfo = TRUE) 
{
  if (missing(def)) {
    stop(paste("def is missing", sep=""))
  }
  if(def=="point"){
    backgType <- "modelBased"
    } 
  if (missing(backgType)) {
    stop(paste("backgType is missing", sep=""))
  }
  if (!(def %in% c("excess", "additional", "relative", "extra", "added", "hybridExc", "hybridAdd", "point"))) {
    stop(paste("Could not recognize def", sep=""))
  }
  if (!(backgType %in% c("modelBased","absolute","hybridSD","hybridPercentile"))) {
    stop(paste("Could not recognize backgType", sep=""))
  }
  
  interval <- match.arg(interval)
  
  # Extract information from model
  EDlist <- object$fct[["edfct"]]  
  parmMat <- object$parmMat
  nCurves <- ncol(parmMat)
  
  # bmrScaledList
  bmrScaledList <- getBmrScaledList(object, bmr, backgType, backg, controlSD, def)
  
  # SINGLE CURVE
  if(nCurves == 1){
    
    # Initialise result matrices
    resMat <- matrix(NA, nrow = 1, ncol = 2, dimnames = list("", c("BMD", "BMDL")))
    bmrScaled <- matrix(bmrScaledList$bmrScaledMat, nrow = 1, ncol = 1, dimnames = list("", "bmrScaled"))
    bmdInterval <- matrix(NA, nrow = 1, ncol = 2, dimnames = list("", c("Lower", "Upper")))
    bmdSE <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("", "SE"))
    
    EDeval <- EDlist(parmMat[, 1], bmrScaled[1,], type = "absolute")
    bmdVal <- EDeval[[1]]
    bmdSEVal <- NA
    
    if(interval == "delta"){
      if(sandwich.vcov){
        varCov <- sandwich(object)
      } else {
        varCov <- vcov(object)
      }
      dBmdVal <- EDeval[[2]] + bmrScaledList$dBmrScaled[,1] / object$fct$derivx(bmdVal, t(parmMat))[1]
      bmdSEVal <- sqrt(dBmdVal %*% varCov %*% dBmdVal)
      intMat <- drc:::confint.basic(matrix(c(bmdVal, bmdSEVal), ncol = 2), 
                                    level = 1-2*(1-level), object$"type", df.residual(object), FALSE)
    } else if(interval == "inv"){
      slope <- drop(ifelse(object$curve[[1]](0)-object$curve[[1]](Inf)>0,"decreasing","increasing"))
      if(is.na(object$curve[[1]](0)-object$curve[[1]](Inf))){
        slope <- drop(ifelse(object$curve[[1]](0.00000001)-object$curve[[1]](100000000)>0,"decreasing","increasing"))
      }
      intMat <- invBmd(object, bmr, level = 1-2*(1-level), slope=slope, backgType=backgType,
                       backg=backg, catLev=NA, extFactor=10, def=def, useSD=useSD, 
                       sandwich.vcov=sandwich.vcov)[,c("BMDL", "BMDU"), drop = FALSE]
    } else if(interval == "profile"){
      tmpVals <- ED(object, bmrScaled, interval = "delta",
                        level = 1-2*(1-level), type = "absolute", vcov. = vcov, display = FALSE)[,c("Lower", "Upper"), drop = FALSE]
      if(backgType %in% c("modelBased", "absolute") & def %in% c("excess", "additional","relative", "extra", "point")
         & substr(object$fct$name,1,2) %in% c("LL", "LN", "W1", "W2")){
        if(missing(profileGridSize)){
          profileGridSize <- 20
        }
        slope <- drop(ifelse(object$curve[[1]](0)-object$curve[[1]](Inf)>0,"decreasing","increasing"))
        if(is.na(object$curve[[1]](0)-object$curve[[1]](Inf))){
          slope <- drop(ifelse(object$curve[[1]](0.00000001)-object$curve[[1]](100000000)>0,"decreasing","increasing"))
        }
        tmpInterval <- bmdProfileCI(object, bmr = bmr, backgType = backgType, def = def, level = level, gridSize = profileGridSize,
                                    bmdEst = bmdVal, lower = tmpVals[,"Lower"], upper = tmpVals[,"Upper"], slope = slope)
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
    if(missing(profileGridSize)){profileGridSize <- 50}
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
        vcMat <- sandwich(object)
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
        
        EDeval <- EDlist(parmChosen, bmrScaled[iCurve,], type = "absolute") 
        bmdVal <- EDeval[[1]]
        dBmdVal <- EDeval[[2]] + bmrScaledList$dBmrScaled[,iCurve] / object$fct$derivx(bmdVal, t(parmMat))[iCurve]
        bmdSEVal <- sqrt(dBmdVal %*% varCov %*% dBmdVal)
        
        if(interval == "delta"){
          intMat <- drc:::confint.basic(matrix(c(bmdVal, bmdSEVal), ncol = 2), 
                                        level = 1-2*(1-level), object$"type", df.residual(object), FALSE)
        }
        
        resMat[iCurve,] <- c(bmdVal, intMat[1,1])
        bmdInterval[iCurve,] <- intMat
        bmdSE[iCurve,] <- bmdSEVal
      }
    } else {
      # CURVES ARE FITTED INDEPENDENTLY
      bmdCall <- function(object){
        bmd(object, bmr, backgType, backg, controlSD,
            def, interval, sandwich.vcov, display = FALSE, level, profileGridSize, profileProgressInfo)
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



