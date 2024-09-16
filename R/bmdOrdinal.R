bmdOrdinal <- function(object, bmr=0.1, backgType = "modelBased", backg = NA, def="excess", interval = "delta", level = 0.95, R = 500, bootType = "nonparametric", display = TRUE, progressInfo = TRUE){
  bmdList <- lapply(object$drmList, function(mod) bmd(mod, bmr = bmr, backgType = backgType, backg = backg, def=def, display=FALSE))
  BMD <- sapply(bmdList, FUN=function(x) x$Results[1])
  
  if(interval %in% c("delta", "profile")){
    # CI <- bmdOrdinalDeltaCI(object = object, bmr = bmr, backgType = backgType, backg = backg, def = def, level = level)
    bmdList <- lapply(object$drmList, function(mod) bmd(mod, bmr = bmr, backgType = backgType, backg = backg, def=def, interval = interval, level = level, display=FALSE))
    CI <- t(sapply(bmdList, FUN=function(x) x$interval))
  } else if(interval == "bootstrap") {
    bootData <- bootDataGenOrdinal(object, R=R, bootType = bootType)
    
    if(progressInfo){
      cat("Performing bootstrap\n")
      pb <- txtProgressBar(min = 0, max = R, style = 3)
    }
    bmdBootVal <- matrix(NA, nrow = R, ncol = length(object$levelsMerged))
    for(i in 1:R){
      modelBoot <- suppressWarnings(try(drmOrdinal(object$levels, object$dose, object$weights, bootData[[i]], object$fct), silent = TRUE))
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
