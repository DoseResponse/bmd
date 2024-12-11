bmdOrdinal <- function(object, bmr=0.1, backgType = "modelBased", backg = NA, def="excess", interval = c("delta", "sandwich", "profile", "bootstrap"), level = 0.95, R = 500, bootType = c("nonparametric", "parametric", "model", "hierarchical"), display = TRUE, progressInfo = TRUE){
  interval <- match.arg(interval)
  if(is.null(object$blocks)){
    bootType <- match.arg(bootType)
  } else if(length(bootType)!=1){
    bootType <- "hierarchical"
  }
  if(!is.null(object$blocks)){
    if(interval %in% c("delta", "profile")){
      cat('interval types "delta" and "profile" not recommended for ordinal dose-response models with blocks. "sandwich" or "bootstrap" with bootType = "hierarchical" is recommended.\n')
    }
    if((interval == "bootstrap") & (bootType != "hierarchical")){
      cat('bootType = "hierarchical" is recommended for ordinal dose-response models with blocks.\n')
    }
  }
  
  if(interval %in% c("delta", "sandwich", "profile")){
    # CI <- bmdOrdinalDeltaCI(object = object, bmr = bmr, backgType = backgType, backg = backg, def = def, level = level)
    bmdList <- lapply(object$drmList, function(mod) bmd(mod, bmr = bmr, backgType = backgType, backg = backg, def=def, interval = interval, level = level, display=FALSE))
    BMD <- sapply(bmdList, FUN=function(x) x$Results[1])
    CI <- t(sapply(bmdList, FUN=function(x) x$interval))
  } else if(interval == "bootstrap") {
    bmdList <- lapply(object$drmList, function(mod) bmd(mod, bmr = bmr, backgType = backgType, backg = backg, def=def, display=FALSE))
    BMD <- sapply(bmdList, FUN=function(x) x$Results[1])
    
    bootData <- bootDataGenOrdinal(object, R=R, bootType = bootType)
    
    if(progressInfo){
      cat("Performing bootstrap\n")
      pb <- txtProgressBar(min = 0, max = R, style = 3)
    }
    bmdBootVal <- matrix(NA, nrow = R, ncol = length(object$levelsMerged))
    for(i in 1:R){
      modelBoot <- suppressWarnings(try(drmOrdinal(levels = object$levels, dose = object$dose, weights = object$weights, data = bootData[[i]], fct = object$fct), silent = TRUE))
      bmdAllBoot <- try(lapply(modelBoot$drmList, function(mod) try(bmd(mod, bmr = bmr, backgType = backgType, backg = backg, def=def, display=FALSE)$Results[,1], silent = TRUE)), silent = TRUE)
      bmdBootVal[i,] <- suppressWarnings(sapply(bmdAllBoot, as.numeric))
      if(progressInfo) setTxtProgressBar(pb, i)
    }
    if(progressInfo) close(pb)
    
    BMDL <- apply(bmdBootVal, 2, quantile, p= c(1-level), na.rm = TRUE)
    BMDU <- apply(bmdBootVal, 2, quantile, p= c(level), na.rm = TRUE)
    # CI <- quantile(bmdBoot, c(1-level, level), na.rm = TRUE)
    CI <- cbind(BMDL, BMDU)
    rownames(CI) <- object$levelsMerged
  } else {
    CI <- NA
  }
  
  resMat <- matrix(data = c(BMD, CI[,1]), ncol = 2, dimnames=list(object$levelsMerged, c("BMD", "BMDL")))
  
  resBMD<-list(Results = resMat,
               interval = CI,
               bmdList = bmdList)
  
  class(resBMD) <- "bmdOrdinal"
  if(display){ print(resMat) }
  invisible(resBMD)
}
