#' Benchmark dose estimation with heterogeneous variance
#' 
#' Estimation of benchmark doses and benchmark dose lower limit based on the
#' hybrid method from dose response model fits with the option to specify a
#' heterogeneous variance structure, where the variance depends on the dose
#' level and/or the fitted values
#' 
#' The aim to provide an R package calculating the benchmark dose (BMD) and the
#' lower limit of the corresponding 95\% confidence interval (BMDL) for
#' continuous and quantal dose-response data for a range of dose-response
#' models based on the available definitions of the benchmark dose concepts.
#' 
#' REFERENCES TO BE ADDED/WRITTEN
#' 
#' @param object dose-response model with a heterogeneous variance structure of
#' class \code{drcHetVar}
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
#' "boot" - BMDL is based on percentile bootstrapping.
#' 
#' "none" - no confidence interval is computed.
#' @param R number of bootstrap samples. Ignored if \code{interval = "none"}
#' @param level numeric value specifying the levle of the confidence interval
#' underlying BMDL. Default is 0.95
#' @param bootType character string specifying the type of bootstrap samples.
#' Options are "nonparametric" (observations are drawn without replacement from
#' the original data set), "semi-parametric" (standardised residuals are drawn
#' with replacement and subsequently rescaled according to the model) and
#' "parametric" (new observations are simulated from the distribution given by
#' the fitted model).
#' @param progressInfo logical. If TRUE, progress info is be printed while
#' bootstrap confidence intervals are estimated. Default is TRUE.
#' @param display logical. If TRUE the results are displayed; otherwise they
#' are not
#' @return A list of five elements: Results contain the estimated BMD and BMDL,
#' bmrScaled is the response value corresponding to the BMD, interval gives the
#' lower (BMDL) and upper (BMDU) end of the confidence interval of BMD
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
#' ryegrass.LL.4.hetVar <- drmHetVar(rootl ~ conc, ~ fitted + I(fitted^2),
#'                                   data = ryegrass, fct = LL.4())
#' plot(ryegrass.LL.4.hetVar)
#' bmdHetVar(ryegrass.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, 
#'           def = "hybridExc", R = 50, level = 0.95, progressInfo = TRUE, 
#'           display = TRUE) # increase R
#' bmdHetVar(ryegrass.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, 
#'           def = "hybridExc", R = 50, level = 0.95, 
#'           bootType = "parametric", progressInfo = TRUE, display = TRUE) # parametric bootstrap
#' 
#' # barley data
#' set.seed(123)
#' barley.LL.4.hetVar <- drmHetVar(weight ~ Dose, ~ fitted + I(fitted^2), data = barley, fct = LL.4())
#' plot(barley.LL.4.hetVar)
#' bmdHetVar(barley.LL.4.hetVar, bmr = 0.1, backgType = "hybridSD", backg = 1, 
#'           def = "hybridExc", R = 50, level = 0.95, progressInfo = TRUE, display = TRUE)
#' 
#' # GiantKelp data
#' set.seed(123)
#' GiantKelp.LL.4.hetVarSq <- drmHetVar(tubeLength ~ dose, ~ fitted + I(fitted^2), 
#'                                      data = GiantKelp, fct = LL.4())
#' plot(GiantKelp.LL.4.hetVarSq)
#' bmdHetVar(GiantKelp.LL.4.hetVarSq, bmr = 0.1, backgType = "hybridSD", backg = 1, 
#'           def = "hybridExc", R = 50, level = 0.95, progressInfo = TRUE, display = TRUE)
#' 
#' GiantKelp.LL.4.hetVarLogSq <- drmHetVar(tubeLength ~ dose, ~ log(dose+1) + I(log(dose+1)^2), 
#'                                         data = GiantKelp, fct = LL.4())
#' plot(GiantKelp.LL.4.hetVarLogSq)
#' bmdHetVar(GiantKelp.LL.4.hetVarLogSq, bmr = 0.1, backgType = "hybridSD", backg = 1, 
#'           def = "hybridExc", R = 50, level = 0.95, progressInfo = TRUE, display = TRUE)
#' 
#' 
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
  
  if(length(bootType) != 1){
    if(!identical(bootType, c("nonparametric", "semiparametric", "parametric"))){
    }
    bootType <- "nonparametric" # default
  }
  if(!bootType %in% c("nonparametric", "semiparametric", "parametric")){
    stop('"bootType" not recognised. Options are: "nonparametric", "semiparametric" and "parametric"')
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
  class(resBMD) <- c("bmdHetVar", "bmd")
  invisible(resBMD) 
}
