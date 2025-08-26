#' Benchmark dose estimation for ordinal dose-response models
#' 
#' Estimation of benchmark doses and benchmark dose lower limit based on model
#' averaging from a user-defined list of ordinal dose-response model fits
#' 
#' The aim to provide an R package calculating the benchmark dose (BMD) and the
#' lower limit of the corresponding 95\% confidence interval (BMDL) for
#' continuous, quantal, quantal and ordinal dose-response data for a range of
#' dose-response models based on the available definitions of the benchmark
#' dose concepts.
#' 
#' Details on the implemented definitions and methods can be found in Crump
#' (2002)
#' 
#' @param modelList A list of ordinal dose-response models of class
#' \code{drcOrdinal}
#' @param modelWeights character string specifying the type of weights used,
#' "AIC" or "BIC", or a numeric vector of the same length as the modelList with
#' user defined weights
#' @param bmr numeric value of benchmark response level for which to calculate
#' the benchmark dose
#' @param backgType character string specifying how the background level is
#' specified. The options are "modelBased" and "absolute".
#' 
#' "modelBased" - the background level is obtained from the model as the level
#' for dose 0: p0 = f(0)
#' 
#' "absolute" - the background level is specified by the user through the backg
#' argument.
#' @param backg numeric value specifying the background level. Defaults to 0
#' for "absolute" background risk.
#' @param def character string specifying the definition of the benchmark dose
#' to use in the calculations. "excess", "additional" and "point" are available
#' for ordinal response.
#' 
#' "excess" - BMR is defined as: BMR = (f(BMD) - p0)/(1 - p0).  Works for
#' binomial response. BMR should be between 0 and 1.
#' 
#' "additional" - BMR is defined as: BMR = f(BMD) - p0.  Works for binomial
#' response. BMR should be between 0 and 1.
#' 
#' "point" - The response level for which to find BMD is directly defined
#' through the BMR level: BMR = f(BMD). Works for binomial, count and
#' continuous response data.
#' @param type specify the model averaging type used for the confidence
#' intervals. Options are "Kang" (confidence intervals computed as weighted
#' averages of individual confidence intervals) and "bootstrap"
#' @param level numeric value specifying the levle of the confidence interval
#' underlying BMDL. Default is 0.95
#' @param R integer specifying the number of data sets resampled from the
#' original data set when computing the confidence interval by bootstrap.
#' @param bootType character string specifying the resampling procedure for the
#' data sets used for bootstrap.
#' 
#' "nonparametric" - Bootstrapping is done by sampling with replacement from
#' the original data set.
#' 
#' "parametric" - Bootstrapping is done by sampling from a multinomial
#' distribution with probabilites given by the number of observations for each
#' level. If all observations for one group are in the same level, shrinkage is
#' used to avoid that the resampling always produces that particular level. In
#' this case data is sampled from a multinomial distribution with probabilities
#' (1/|K|^2)/(N_i + 1/|K|) for the levels with 0 observations, and (N_i +
#' 1/|K|^2)/(N_i + 1/|K|), where N_i is the number of observatiions for the
#' particular group, and |K| is the number of levels.
#' 
#' All resampling is done within the dose values.
#' @param display logical. If TRUE the results are displayed; otherwise they
#' are not
#' @param progressInfo logical. If TRUE, a progressbar is displayed when
#' calculating bootstrap confidence intervals
#' @return A list of four elements: Results contain the estimated BMD and BMDL,
#' interval gives the lower (BMDL) and upper (BMDU) end of the confidence
#' interval of BMD.
#' @author Signe M. Jensen and Jens Riis Baalkilde
#' @references Budtz-Jorgensen, E., Keiding, N., and Grandjean, P. (2001)
#' Benchmark Dose Calculation from Epidemiological Data, \emph{Biometrics}
#' \bold{57}, 698--706.
#' 
#' Crump, K. (2002) Critical Issues in Benchmark Calculations from Continuous
#' Data, \emph{Critical Reviews in Toxicology} \bold{32}, 133--153.
#' @keywords models nonlinear
#' @examples
#' 
#' library(drc)
#' library(drcData)
#' data(guthion)
#' 
#' guthionS <- subset(guthion, trt == "S")
#' 
#' guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), 
#'                           weights = "total", dose = "dose", data = guthionS, fct = LL.2())
#' guthionS.LN <- drmOrdinal(levels = c("alive", "moribund", "dead"), 
#'                           weights = "total", dose = "dose", data = guthionS, fct = LN.2())
#' guthionS.W1 <- drmOrdinal(levels = c("alive", "moribund", "dead"),
#'                           weights = "total", dose = "dose", data = guthionS, fct = W1.2())
#' guthionS.W2 <- drmOrdinal(levels = c("alive", "moribund", "dead"),
#'                           weights = "total", dose = "dose", data = guthionS, fct = W2.2())
#' 
#' bmdOrdinalMA(list(guthionS.LL, guthionS.LN, guthionS.W1, guthionS.W2), 
#'              modelWeights = "AIC", bmr = 0.1, 
#'              backgType = "modelBased", def = "excess", type = "Kang")
#' bmdOrdinalMA(list(guthionS.LL, guthionS.LN, guthionS.W1, guthionS.W2), 
#'              modelWeights = "AIC", bmr = 0.1, 
#'              backgType = "modelBased", def = "excess", type = "bootstrap", R = 50)
#' @export
bmdOrdinalMA <- function(modelList, modelWeights = c("AIC", "BIC"), bmr, backgType = c("modelBased", "absolute"), backg = NA, def = c("excess", "additional", "point"), type = c("bootstrap", "Kang"), level = 0.95, R = 500, bootType = c("nonparametric", "parametric"), display = TRUE, progressInfo = TRUE){
  # assertions
  if(!all(sapply(modelList, function(object) inherits(object, "drcOrdinal")))){
    stop('"modelList" must be a list of ordinal dose-response models of type "drcOrdinal"')
  }
  
  bootType <- match.arg(bootType)
  
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
