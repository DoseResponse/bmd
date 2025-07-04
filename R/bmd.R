#' Benchmark dose estimation
#' 
#' Estimation of benchmark doses and benchmark dose lower limit from dose
#' response model fits
#' 
#' The aim to provide an R package calculating the benchmark dose (BMD) and the
#' lower limit of the corresponding 95\% confidence interval (BMDL) for
#' continuous and quantal dose-response data for a range of dose-response
#' models based on the available definitions of the benchmark dose concepts.
#' 
#' Details on the implemented definitions and methods can be found in Crump
#' (2002)
#' 
#' @param object object of class \code{drc}
#' @param bmr numeric value of benchmark response level for which to calculate
#' the benchmark dose
#' @param backgType character string specifying how the background level is
#' specified. For binomial data the options are "modelBased" and "absolute".
#' For continuous data the options are "modelBased","absolute", "hybridSD" and
#' "hybridPercentile". For count data (Poisson, negbin1 or negbin2) the options
#' are "modelBased" and "absolute".
#' 
#' "modelBased" - the background level is obtained from the model as the level
#' for dose 0: p0 = f(0)
#' 
#' "absolute" - the background level is specified by the user through the backg
#' argument: p0 = backg for binomial response and for the "relative", "extra"
#' and "added" definition for continuous or count response data.  p0 = 1 -
#' phi((back - f(0))/sigma) for "hybridExc" and "hybridAdd" definitions.
#' 
#' "hybridSD" - the background risk is specified by the user in terms of number
#' of SDs from the mean of the control group.  p0 = 1 - phi(((backg*sigma +
#' f(0)) - f(0))/sigma) = 1 - phi(backg), where phi is the normal distribution
#' function and sigma is the SD for the control group.
#' 
#' "hybridPercentile" - the background risk is specified by the user in terms
#' of percentile from the control group distribution (assuming a normal
#' distribution).  p0 = 1 - phi((x0 - f(0))/sigma) = 1 - backg.  where x0 is
#' the level for which the response is considered adverse, phi is the normal
#' distribution function and sigma is the SD for the control group
#' 
#' If not specified, it will be set to "modelBased" for all def excluding
#' "hybridExc" and "hybridAdd", for these the default is "hybridSD".
#' @param backg numeric value specifying the background level. Defaults to 0
#' for "absolute" background risk for binomial response (1 for decreasing
#' dose-response models), 2 SD for "hybridSD" background and 0.9 for
#' "hybridPercentile"
#' @param controlSD numeric value specifying the standard deviation of the
#' control group (used in the hybrid approach). If not specified the SD for the
#' control group will be estimated as the square root of the residual variance
#' (assuming variance homogeneity).
#' @param def character string specifying the definition of the benchmark dose
#' to use in the calculations. "excess", "additional" and "point" are for
#' binomial response. "relative", "extra", "added", "hybridExc" (excess
#' hybrid), "hybridAdd" (additional hybrid), and "point" are for continuous
#' response. "relative", "extra", and "point" are for count response data.
#' 
#' "excess" - BMR is defined as: BMR = (f(BMD) - p0)/(1 - p0).  Works for
#' binomial response. BMR should be between 0 and 1.
#' 
#' "additional" - BMR is defined as: BMR = f(BMD) - p0.  Works for binomial
#' response. BMR should be between 0 and 1.
#' 
#' "point" - The response level for which to find BMD is directly defined
#' through the BMR level: BMR = f(BMD). Works for binomial, count and
#' continuous response data
#' 
#' "relative" - BMR is defined as: BMR = (f(BMD) - p0)/p0.  Works for
#' continuous and count response data
#' 
#' "extra" - BMR is defined as: BMR = (f(BMD) - p0)/(f(Inf) - p0).  Works for
#' continuous and count response data
#' 
#' "added" - BMR is defined as: BMR= f(BMD) + p0.  Works for continuous
#' response
#' 
#' "hybridExc" - BMR is defined as: BMR = (1 - phi((x0 - f(BMD))/sigma) - p0)/
#' (1- p0), where x0 is the level for which the response is considered adverse,
#' phi is the normal distribution function and sigma is the SD for the control
#' group.  Works for continuous response
#' 
#' "hybridAdd" - BMR is defined as: BMR = 1 - phi((x0 - f(BMD))/sigma) - p0,
#' where x0 is the level for which the response is considered adverse, phi is
#' the normal distribution function and sigma is the SD for the control group.
#' Works for continuous response
#' @param respTrans if the dose-response model is fitted with a transformed
#' response, specifying this option ensures that the background level for the
#' BMD is computed on the original scale.  Options include "none" (default),
#' "log" (natural logarithm) and "sqrt"(square root).
#' @param interval character string specifying the type of confidence interval
#' to use: "delta" (default), "inv", "profile" or "profileGrid"
#' 
#' "delta" - BMDL is based on the lower limit of a Wald confidence interval
#' based on the delta method
#' 
#' "inv" - BMDL is based on inverse regression, that is the dose associated
#' with the upper limit of the pointwise confidence interval of the
#' dose-response curve. Currently not available for a transformed response
#' variable.
#' 
#' "profile" - A profile confidence interval is computed by reparametrising the
#' model. Only available for \code{backgType} in ("modelBased", "absolute"),
#' \code{def} in ("excess", "additional","relative", "extra", "point") and
#' dose-response models with model functions of types Log-Logistic, Log-Normal
#' and Weibull 1 and 2. The limits of the confidence interval are searched for
#' within a grid of size specified by the \code{profileGridSize} argument,
#' however the values of the confidence interval are not limited to the values
#' on the grid.
#' 
#' "profileGrid" - A profile confidence interval is computed by creating a
#' multi-dimensional grid of parameter values for the model. Then, a confidence
#' region for the parameters in the model is estimated, and the BMD is
#' estimated for all combinations of parameter values within the region. The
#' interval for the BMD is then computed as the minimum and maximum values of
#' BMD values within the region. Note that for models with several parameters
#' (3-5), is might take a while to compute the confidence interval using this
#' method. The number of grid points for each parameter can be specified by the
#' \code{profileGridSize} argument.
#' @param sandwich.vcov logical. If TRUE BMDL is based on confidence intervals
#' based on robust standard errors from the sandwich covariance matrix
#' @param display logical. If TRUE the results are displayed; otherwise they
#' are not
#' @param level numeric value specifying the levle of the confidence interval
#' underlying BMDL. Default is 0.95
#' @param profileGridSize integer specifying the number of grid points used for
#' each parameter, when \code{interval} is either "profile" or
#' "profileGridSearch". Defaults to 20 if not specified.
#' @param profileProgressInfo logical. If TRUE, progress info is be printed
#' while computed profile confidence intervals using the \code{profileGrid}
#' method. Default is TRUE.
#' @return A list of four elements: Results contain the estimated BMD and BMDL,
#' bmrScaled is the response value corresponding to the BMD, interval gives the
#' lower (BMDL) and upper (BMDU) end of the confidence interval of BMD, SE
#' gives the standard error of BMD.
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
#' 
#' ## Fitting log-logistic two-parameter model to binomial data
#' deguelin.m1 <- drm(r/n~dose, weights=n, data=deguelin, fct=LL.2(), type="binomial")
#' 
#' ## BMD for 5% additional risk with estimated background risk
#' bmd(deguelin.m1, 0.05, backgType = "modelBased", def = "additional")
#' 
#' ## BMD for 10% additional risk with 2% background risk
#' bmd(deguelin.m1, 0.1, backg = 0.02 , backgType = "absolute", def = "additional")
#' 
#' ## BMD for 5% excess risk and background 0
#' bmd(deguelin.m1, 0.05, backg = 0, backgType = "absolute", def = "excess")
#' 
#' ## Dose resulting in 12% risk
#' bmd(deguelin.m1, 0.12, def = "point")
#' 
#' ## Benchmark doses for a continuous response
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#' 
#' ## BMD as the dose resulting in a 5% change relative to the mean background level
#' bmd(ryegrass.m1, 0.05, backgType = "modelBased", def = "relative", display = TRUE)
#' 
#' ## BMD using the hybrid method, background risk is 2 SD, hybrid definition using excess risk
#' bmd(ryegrass.m1, 0.05, backg = 2, backgType = "hybridSD", def = "hybridAdd", display = TRUE)
#' 
#' 
#' @export
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
        intMat <- drc:::confint.basic(matrix(c(bmdVal, bmdSEVal), ncol = 2), 
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
          intMat <- drc:::confint.basic(matrix(c(bmdVal, bmdSEVal), ncol = 2), 
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



