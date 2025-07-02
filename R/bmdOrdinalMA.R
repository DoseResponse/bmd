bmdOrdinalMA <- function(modelList, modelWeights = c("AIC", "BIC"), bmr=0.1, backgType = "modelBased", backg = NA, def="excess", type = c("bootstrap", "Kang"), level = 0.95, R = 500, bootType = "nonparametric", display = TRUE, progressInfo = TRUE){
  # assertions
  if(!all(sapply(modelList, function(object) inherits(object, "drcOrdinal")))){
    stop('"modelList" must be a list of ordinal dose-response models of type "drcOrdinal"')
  }
  
  # bmdEstimates on all models
  bmdList <- lapply(modelList, function(object) bmdOrdinal(object, bmr=bmr, backgType = backgType, backg = backg, def=def, interval = "delta", display = FALSE))
  
  # modelWeights
  if(missing(modelWeights)){
    stop("missing argument \"modelWeights\". Options are \"AIC\", \"BIC\" or a numeric vector of same lenght as modelList.")
  } else if(length(modelWeights) > 1){
    if(!identical(length(modelWeights), length(modelList))){ 
      stop("misspecified argument \"modelWeights\". Options are \"AIC\", \"BIC\" or a numeric vector of same lenght as modelList.")
    } else {
      modelWeights0 <- modelWeights / sum(modelWeights)
    }
  } else if(length(modelWeights) == 1){
    if(identical(modelWeights, "AIC")){
      AICVals <- sapply(modelList, AIC)
      modelWeights0 <- exp(-1/2*(AICVals - min(AICVals)))/sum(exp(-1/2*(AICVals - min(AICVals))))
    } else if(identical(modelWeights, "BIC")){
      BICVals <- sapply(modelList, BIC)
      modelWeights0 <- exp(-1/2*(BICVals - min(BICVals)))/sum(exp(-1/2*(BICVals - min(BICVals))))
    }
  }
  
  if(all(type == c("bootstrap","Kang")) | sum(type[1] == c("bootstrap","Kang")) != 1){
    stop('Specify model averaging type. Options are "bootstrap" and "Kang"')
  }
  
  if(identical(type,"Kang")){
    maBMD <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$Results[,1])))
    maBMDL <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$interval[,1])))
    maBMDU <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$interval[,2])))
  } else if(identical(type, "bootstrap")){
    maBMD <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$Results[,1])))
    
    # bootstrap
    bootData <- bootDataGenOrdinal(modelList[[1]],R=R, bootType = bootType)
    
    bmdMAboot <- function(data){
      bootModelList <- lapply(modelList, function(model) try(
        eval(substitute(drmOrdinal(levels = levels0, dose = dose0, weights = weights0, blocks = blocks0, data = bootData[[i]], fct = model$fct),
                        list(levels0 = modelList[[1]]$levels,
                             dose0 = modelList[[1]]$dose,
                             weights0 = modelList[[1]]$weights,
                             blocks0 = modelList[[1]]$blocks))),
        silent = TRUE))
      
      modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
      
      bootModelList <- bootModelList[!modelConvergenceError]
      
      if(length(modelWeights) > 1){
        bootModelWeights0 <- modelWeights0[!modelConvergenceError]
      } else if(length(modelWeights) == 1){
        if(identical(modelWeights, "AIC")){
          BootAICVals <- sapply(bootModelList, AIC)
          bootModelWeights0 <- exp(-1/2*(BootAICVals - min(BootAICVals)))/sum(exp(-1/2*(BootAICVals - min(BootAICVals))))
        } else if(identical(modelWeights, "BIC")){
          BootBICVals <- sapply(bootModelList, BIC)
          bootModelWeights0 <- exp(-1/2*(BootBICVals - min(BootBICVals)))/sum(exp(-1/2*(BootBICVals - min(BootBICVals))))
        }
      }
      
      bootBmdList <- lapply(bootModelList, function(object) bmdOrdinal(object, bmr=bmr, backgType = backgType, def=def, interval = "delta", display = FALSE))
      
      bootMaBMD <- colSums(bootModelWeights0 * t(sapply(bootBmdList, function(x) x$Results[,1])))
      bootMaBMD
    }
    
    bootBmdEst <- matrix(NA, nrow = length(bootData), ncol = length(modelList[[1]]$levelsMerged))
    
    if(progressInfo){
      cat("Performing bootstrap\n")
      pb <- txtProgressBar(min = 0, max = R, style = 3)
    }
    
    for(i in 1:length(bootData)){
      bootBmdEst[i,] <- bmdMAboot(bootData[[i]])
      if(progressInfo) setTxtProgressBar(pb, i)
    }
    if(progressInfo) close(pb)
    
    bootBmdError <- apply(bootBmdEst, 1, function(x) any(is.na(x)))
    boot0 <- bootBmdEst[!bootBmdError,]
    
    if(length(boot0) == 0){ 
      maBMDL <- NA 
      maBMDU <- NA
    } else {
      maBMDL <- apply(boot0, 2, quantile, p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
      maBMDU <- apply(boot0, 2, quantile, p=c(level), na.rm = TRUE)
    }
    # end bootstrap
  }
  
  resMat<-matrix(c(maBMD,maBMDL), nrow = length(modelList[[1]]$levelsMerged), ncol = 2, byrow = FALSE)
  colnames(resMat) <- c("BMD_MA", "BMDL_MA")
  rownames(resMat) <- modelList[[1]]$levelsMerged
  
  intMat<-matrix(c(maBMDL,maBMDU), nrow = length(modelList[[1]]$levelsMerged), ncol = 2, byrow = FALSE)
  colnames(intMat) <- c("BMDL_MA", "BMDU_MA")
  rownames(intMat) <- modelList[[1]]$levelsMerged
  
  used.Boot<-ifelse(identical(type,"bootstrap")|identical(type,"Bootstrap"),
                    floor(length(boot0)/length(modelList[[1]]$levelsMerged)),NA)
  
  resBMD<-list(Results = resMat,
               Boot.samples.used = used.Boot,
               interval = intMat,
               modelWeights = modelWeights0,
               bmdList = bmdList)
  
  if(display){
    print(resMat)
  }
  
  class(resBMD) <- "bmdOrdinal"
  invisible(resBMD)
}

# old bmdOrdinalMA function -----------------------------------------------

# bmdOrdinalMA <- function(modelList, bmr=0.1, backgType = "modelBased", def="excess", level = 0.95, R = 500, bootType = "nonparametric", MAType = "postBmdEst", display = TRUE){
#   MABmdEst <- function(modelList){
#     bmdAllMods <- lapply(modelList, function(object) bmdOrdinal(object, bmr=bmr, backgType = backgType, def=def, interval = "none", display = FALSE))
#     
#     AICVals <- sapply(modelList, AIC)
#     AICWeights <- exp(-1/2*(AICVals - min(AICVals)))/sum(exp(-1/2*(AICVals - min(AICVals))))
#     
#     if(MAType == "postBmdEst"){
#       BMD.tmp <- sapply(bmdAllMods, function(object) object$Results[1])
#       bmdVal <- sum(BMD.tmp*AICWeights)
#     } else if (MAType == "preBmdEst"){
#       bmdEachCat <- sapply(1:length(modelList[[1]]$drmList), function(cat_i){
#         bmdAllMods0 <- sapply(bmdAllMods, 
#                               function(object){
#                                 object$bmdList[[cat_i]]$Results[1]
#                               })
#         sum(AICWeights*bmdAllMods0)
#       })
#       
#       bmdVal <- mean(bmdEachCat)
#     }
#     bmdVal
#   }
#   
#   BMD <- MABmdEst(modelList)
#   
#   # bootstrap
#   bootData <- bootDataGenOrdinal(modelList[[1]],R=R, bootType = bootType)
#   
#   bmdBoot <- numeric(R)
#   for(i in 1:R){
#     modelListBoot <- lapply(modelList, 
#                             function(object) suppressWarnings(try(drmOrdinal(object$levels, object$dose, object$weights, bootData[[i]], object$fct), silent = TRUE))
#                             )
#     bmdBoot[i] <- suppressWarnings(as.numeric(try(MABmdEst(modelListBoot), silent = TRUE)))
#   }
#   CI <- quantile(bmdBoot, c(1-level, level), na.rm = TRUE)
#   Results = c(BMD = BMD, BMDL = CI[1])
#   
#   bmdRes <- list(Results = Results, 
#                  interval = CI
#                  )
#   if(display) print(Results)
#   invisible(bmdRes)
# }
