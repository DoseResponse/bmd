#' Benchmark dose estimation with heterogeneous variance based on model
#' averaging (MA)
#' 
#' Estimation of benchmark doses and benchmark dose lower limit based on the
#' hybrid method from a list of dose response model fits with the option to
#' specify a heterogeneous variance structure, where the variance depends on
#' the dose level and/or the fitted values
#' 
#' The aim to provide an R package calculating the benchmark dose (BMD) and the
#' lower limit of the corresponding 95\% confidence interval (BMDL) for
#' continuous and quantal dose-response data for a range of dose-response
#' models based on the available definitions of the benchmark dose concepts.
#' 
#' REFERENCES TO BE ADDED/WRITTEN
#' 
#' @param modelList a list of dose-response models with a heterogeneous
#' variance structure of class \code{drcHetVar}
#' @param modelWeights character string specifying the type of weights used,
#' "AIC" or "BIC", or a numeric vector of the same length as the modelList with
#' user defined weights
#' @param bmr numeric value of benchmark response level for which to calculate
#' the benchmark dose
#' @param backgType character string specifying how the background level is
#' specified. The options are "absolute", "hybridSD" and "hybridPercentile".
#' 
#' "absolute" - the background level is specified by the user through the backg
#' argument: p0 = 1 - phi((back - f(0))/sigma(0)) for "hybridExc" and
#' "hybridAdd" definitions.
#' 
#' "hybridSD" - the background risk is specified by the user in terms of number
#' of SDs from the mean of the control group.  p0 = 1 - phi(((backg*sigma(0) +
#' f(0)) - f(0))/sigma(0)) = 1 - phi(backg), where phi is the normal
#' distribution function and sigma(0) is the SD for the control group.
#' 
#' "hybridPercentile" - the background risk is specified by the user in terms
#' of percentile from the control group distribution (assuming a normal
#' distribution).  p0 = 1 - phi((x0 - f(0))/sigma(0)) = 1 - backg.  where x0 is
#' the level for which the response is considered adverse, phi is the normal
#' distribution function and sigma(0) is the SD for the control group
#' @param backg numeric value specifying the background level. Defaults to 2 SD
#' for "hybridSD" background and 0.9 for "hybridPercentile"
#' @param def character string specifying the definition of the benchmark dose
#' to use in the calculations. "hybridExc" (excess hybrid), "hybridAdd"
#' (additional hybrid), available.
#' 
#' "hybridExc" - BMR is defined as: BMR = (1 - phi((x0 - f(BMD))/sigma(BMD)) -
#' p0)/ (1- p0), where x0 is the level for which the response is considered
#' adverse, phi is the normal distribution function and sigma(BMD) is the SD at
#' the benchmark dose.
#' 
#' "hybridAdd" - BMR is defined as: BMR = 1 - phi((x0 - f(BMD))/sigma(BMD)) -
#' p0, where x0 is the level for which the response is considered adverse, phi
#' is the normal distribution function and sigma(BMD) is the SD at the
#' benchmark dose.
#' @param interval character string specifying the type of confidence interval
#' to use: "boot" (default) or "none"
#' 
#' "boot" - BMDL is based on nonparametric percentile bootstrapping.
#' 
#' "none" - no confidence interval is computed.
#' @param R number of bootstrap samples. Ignored if \code{interval = "none"}
#' @param level numeric value specifying the levle of the confidence interval
#' underlying BMDL. Default is 0.95
#' @param progressInfo logical. If TRUE, progress info is be printed while
#' bootstrap confidence intervals are estimated. Default is TRUE.
#' @param display logical. If TRUE the results are displayed; otherwise they
#' are not
#' @return A list of four elements: Results contain the estimated BMD and BMDL,
#' Boot.samples.used gives the number of boot samples that resulted in
#' succesful estimations and were accordingly used in the estimation of BMDL
#' (and BMDU), Interval gives BMDL and BMDU, which is identical to the
#' confidence interval for the percentile interval approach, and modelWeights
#' includes the estimated weights.
#' @author Signe M. Jensen and Jens Riis Baalkilde
#' @keywords models nonlinear
#' @examples
#' 
#' library(drc)
#' library(drcData)
#' library(bmd)
#' # install.packages("gridExtra") # OPTIONAL - USED FOR PLOTTING A drcHetVar OBJECT.
#' 
#' # ryegrass data
#' set.seed(123)
#' ryegrass.hetVar.list <- list(
#'   drmHetVar(rootl ~ conc, ~ fitted + I(fitted^2), data = ryegrass, fct = LL.4()),
#'   drmHetVar(rootl ~ conc, ~ fitted + I(fitted^2), data = ryegrass, fct = LN.4()),
#'   drmHetVar(rootl ~ conc, ~ fitted + I(fitted^2), data = ryegrass, fct = W1.4()),
#'   drmHetVar(rootl ~ conc, ~ fitted + I(fitted^2), data = ryegrass, fct = W2.4()))
#' bmdHetVarMA(ryegrass.hetVar.list, modelWeights = "AIC", bmr = 0.1, backgType = "hybridPercentile",
#'             backg = 0.1, def = "hybridExc", R = 100, level = 0.95)
#' bmdHetVarMA(ryegrass.hetVar.list, modelWeights = c(0.4, 0.2, 0.1, 0.3), bmr = 0.1, 
#'             backgType = "hybridPercentile", backg = 0.1, 
#'             def = "hybridExc", R = 50, level = 0.95) # user-defined weights
#' 
#' # barley data
#' set.seed(123)
#' barley.hetVar.list <- list(drmHetVar(weight ~ Dose, ~ fitted + I(fitted^2), 
#'                                      data = barley, fct = LL.4()),
#'                            drmHetVar(weight ~ Dose, ~ fitted + I(fitted^2), 
#'                                      data = barley, fct = W2.4()))
#' bmdHetVarMA(barley.hetVar.list, modelWeights = "AIC", bmr = 0.1, backgType = "hybridSD", backg = 2,
#'             def = "hybridExc", R = 50, level = 0.95, progressInfo = TRUE, display = TRUE)
#' 
#' # GiantKelp data
#' set.seed(123)
#' GiantKelp.hetVar.list <- list(
#'   drmHetVar(tubeLength ~ dose, ~ fitted + I(fitted^2), data = GiantKelp, fct = LL.4()),
#'   drmHetVar(tubeLength ~ dose, ~ log(dose+1) + I(log(dose+1)^2), data = GiantKelp, fct = LL.4()))
#' bmdHetVarMA(GiantKelp.hetVar.list, modelWeights = "AIC", bmr = 0.1, backgType = "hybridSD",
#'             backg = 1, def = "hybridExc", R = 50, level = 0.95, progressInfo = TRUE, 
#'             display = TRUE)
#' 
#' 
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
