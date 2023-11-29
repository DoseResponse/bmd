bmdOrdinal <- function(object, bmr=0.1, backgType = "modelBased", backg = NA, def="excess", interval = "delta", level = 0.95, R = 500, bootType = "nonparametric", display = TRUE){
  bmdList <- lapply(object$drmList, function(mod) bmd(mod, bmr = bmr, backgType = backgType, backg = backg, def=def, display=FALSE))
  BMD <- mean(sapply(bmdList, FUN=function(x) x$Results[1]))
  
  if(interval == "delta"){
    CI <- bmdOrdinalDeltaCI(object = object, bmr = bmr, backgType = backgType, backg = backg, def = def, level = level)
  } else if(interval == "bootstrap") {
    bootData <- bootDataGenOrdinal(object, R=R, bootType = bootType)
    
    bmdBoot <- numeric(R)
    for(i in 1:R){
      modelBoot <- suppressWarnings(try(drmOrdinal(object$levels, object$dose, object$weights, bootData[[i]], object$fct), silent = TRUE))
      bmdAllBoot <- lapply(modelBoot$drmList, function(mod) try(bmd(mod, bmr = bmr, backgType = backgType, backg = backg, def=def, display=FALSE)$Results[1], silent = TRUE))
      bmdBoot[i] <- mean(as.numeric(bmdAllBoot))
    }
    CI <- quantile(bmdBoot, c(1-level, level), na.rm = TRUE)
  } else {
    CI <- NA
  }
  
  resMat <- matrix(data = c(BMD, CI[1]), ncol = 2, dimnames=list(NULL, c("BMD", "BMDL")))
  
  resBMD<-list(Results = resMat,
               interval = CI,
               bmdList = bmdList)
  
  class(resBMD) <- "bmdOrdinal"
  if(display){ print(resMat) }
  invisible(resBMD)
}
