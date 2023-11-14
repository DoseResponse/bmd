bmdOrdinalMA <- function(modelList, bmr=0.1, backgType = "modelBased", def="excess", level = 0.95, R = 500, bootType = "nonparametric", MAType = "postBmdEst", display = TRUE){
  MABmdEst <- function(modelList){
    bmdAllMods <- lapply(modelList, function(object) bmdOrdinal(object, bmr=bmr, backgType = backgType, def=def, interval = "none", display = FALSE))
    
    AICVals <- sapply(modelList, AIC)
    AICWeights <- exp(-1/2*(AICVals - min(AICVals)))/sum(exp(-1/2*(AICVals - min(AICVals))))
    
    if(MAType == "postBmdEst"){
      BMD.tmp <- sapply(bmdAllMods, function(object) object$Results[1])
      bmdVal <- sum(BMD.tmp*AICWeights)
    } else if (MAType == "preBmdEst"){
      bmdEachCat <- sapply(1:length(modelList[[1]]$drmList), function(cat_i){
        bmdAllMods0 <- sapply(bmdAllMods, 
                              function(object){
                                object$bmdList[[cat_i]]$Results[1]
                              })
        sum(AICWeights*bmdAllMods0)
      })
      
      bmdVal <- mean(bmdEachCat)
    }
    bmdVal
  }
  
  BMD <- MABmdEst(modelList)
  
  # bootstrap
  bootData <- bootDataGenOrdinal(modelList[[1]],R=R, bootType = bootType)
  
  bmdBoot <- numeric(R)
  for(i in 1:R){
    modelListBoot <- lapply(modelList, 
                            function(object) suppressWarnings(try(drmOrdinal(object$categories, object$dose, object$weights, bootData[[i]], object$fct), silent = TRUE))
                            )
    bmdBoot[i] <- suppressWarnings(as.numeric(try(MABmdEst(modelListBoot), silent = TRUE)))
  }
  CI <- quantile(bmdBoot, c(1-level, level), na.rm = TRUE)
  Results = c(BMD = BMD, BMDL = CI[1])
  
  bmdRes <- list(Results = Results, 
                 interval = CI
                 )
  if(display) print(Results)
  invisible(bmdRes)
}
