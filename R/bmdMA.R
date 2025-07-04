#' Model-averaged benchmark dose estimation
#' 
#' Estimation of benchmark doses and benchmark dose lower limit based on model
#' averaging from a user-defined list of dose response model fits
#' 
#' This package project is still under development. The aim to provide an R
#' package calculating the benchmark dose (BMD) and the lower limit of the
#' corresponding 95\% confidence interval (BMDL) for continuous and quantal
#' dose-response data for a range of dose-response model based on the available
#' definitions of the benchmark dose concepts.
#' 
#' Details on the implemented definitions and methods can be found in Crump
#' (2002)
#' 
#' Bootstrapping with the argument boot = "nonparametric" is done by sampling
#' with replacement from the original data set. Bootstrapping with the argument
#' boot = "parametric" is done by sampling from norm(mean(Y_i),sd(Y_0)),
#' assuming equal variance between groups, in case of continuous data. For
#' binomial data, each bootstrap data set is sampled from binom(N_i,Y_i/N_i).
#' In case of Y_i = 0 or Y_i = N_i shrinkage is used to avoid that the
#' resampling always produces 0 or 1, respectively. In this case data is
#' sampled from binom(N_i,(Y_i+1/3)/(N_i+1/3)).
#' 
#' All sampling is made within dose groups.
#' 
#' @param modelList list of models of class \code{drc}
#' @param modelWeights character string specifying the type of weights used,
#' "AIC", "BIC" or "Stack" (Baalkilde, J. R., Hansen, N. R., and Jensen, S. M.,
#' 2025), or a vector of the same length as the modelList with user defined
#' weights
#' @param bmr numeric value of benchmark response level for which to calculate
#' the benchmark dose
#' @param backgType character string specifying how the background level is
#' specified. For binomial data the options are "modelBased" and "absolute".
#' For continuous data the options are "absolute", "hybridSD" and
#' "hybridPercentile"
#' 
#' "modelBased" - the background risk is obtained from the model as the risk
#' for dose 0: p0 = f(0)
#' 
#' "absolute" - the background risk is specified by the user through the backg
#' argument: p0 = backg for binomial response and for the "relative", "extra"
#' and "added" definition for continuous response.  p0 = 1 - phi((back -
#' f(0))/sigma) for "hybridExc" and "hybridAdd" definitions.
#' 
#' "hybridSD" - the background risk is specified by the user in terms of number
#' of SDs from the mean of the control group.  p0 = 1 - phi(((backg*sigma +
#' f(0)) - f(0))/sigma) = 1 - phi(backg), where phi is the normal distribution
#' function and sigma is the SD for the control group.  "hybridPercentile" -
#' the background risk is specified by the user in terms of percentile from the
#' control group distribution (assuming a normal distribution).  p0 = 1 -
#' phi((x0 - f(0))/sigma) = 1 - backg.  where x0 is the level for which the
#' response is considered adverse, phi is the normal distribution function and
#' sigma is the SD for the control group
#' @param backg numeric value specifying the background level. Defaults to 0
#' for "absolute" background risk for binomial response (1 for decreasing
#' dose-response models), 2 SD for "hybridSD" background and 0.9 for
#' "hybridpercentile"
#' @param def character string specifying the definition of the benchmark dose
#' to use in the calculations. "excess" , "additional" and "point" are for
#' binomial response whereas "relative", "extra", "added", "hybridExc" (excess
#' hybrid), "hybridAdd" (additional hybrid), and "point" are for continuous
#' response
#' 
#' "excess" - BMR is defined as: BMR = (f(BMD) - p0)/(1 - p0).  Works for
#' binomial response. BMR should be between 0 and 1.
#' 
#' "additional" - BMR is defined as: BMR = f(BMD) - p0.  Works for binomial
#' response. BMR should be between 0 and 1.
#' 
#' "point" - The response level for which to find BMD is directly defined
#' through the BMR level: BMR = f(BMD). Works for both binomial and continuous
#' response
#' 
#' "relative" - BMR is defined as: BMR = (f(BMD) - p0)/p0.  Works for
#' continuous response
#' 
#' "extra" - BMR is defined as: BMR = (f(BMD) - p0)/(f(Inf) - p0).  Works for
#' continuous response
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
#' to use for the individual models (no effect when type = "bootstrap" or
#' "curve"): "delta" (default), "sandwich", "inv" or "profile"
#' 
#' "delta" - BMDL is based on the lower limit of a Wald confidence interval
#' based on the delta method
#' 
#' "sandwich" - BMDL is based on the lower limit of a Wald confidence interval
#' based on the delta method where the sandwich covariance matrix is used
#' 
#' "inv" - BMDL is based on inverse regression, that is the dose associated
#' with the upper limit of the point wise confidence interval of the
#' dose-response curve. Not possible if type = "Buckland"
#' 
#' "profile" - BMDL is based on the lower limit of a profile likelihood
#' confidence interval. Not possible if type = "Buckland"
#' @param type character string specifying how to estimate BMD and BMDL:
#' "curve", "bootstrap", "Kang" or "Buckland"
#' 
#' "curve" -
#' 
#' "bootstrap" -
#' 
#' "Kang" -
#' 
#' "Buckland" -
#' @param bootstrapType character string indicating type of bootstrap sampling
#' to be used if type="bootstrap" or "curve". "nonparametric" (default), or
#' "parametric" (see details)
#' @param bootInterval character string indicating how to estimate the
#' bootstrap confidence intervals used to find BMDL type="bootstrap".
#' "percentile" (default) or "BCa" (Bias corrected and adjusted)
#' @param R number of bootstrap samples to use. Default is 1000
#' @param level numeric value specifying the level of the confidence interval
#' underlying BMDL. Default is 0.95. Note that the full CI is also computed.
#' The full CI has level 1-2*(1-level)
#' @param stackingSeed integer or NULL: Random seed to use in the data split in
#' the estimation of the Stacking Weights when \code{modelWeights = "Stack"}.
#' The global seed is reset to the initial value after estimation of the
#' weights, so this option does not interfere with a globally set seed.
#' @param stackingSplits integer or "LOO": When \code{modelWeights = "Stack"},
#' the Stacking Weights are estimated, which are based on V-fold
#' cross-validation. The stackingSplits argument sets the number V of data
#' splits used in the cross validation. The "LOO" (Leave one out) is a shortcut
#' to setting V equal to the number of observations in the data set.
#' @param display logical. If TRUE the results are displayed; otherwise they
#' are not
#' @param progressInfo logical. If TRUE, a progressbar is shown when computing
#' bootstrap confidence intervals.
#' @return A list of five elements: Results contain the estimated BMD and BMDL,
#' Boot.samples.used gives the number of boot samples that resulted in
#' succesful estimations and were accordingly used in the estimation of BMDL
#' (and BMDU), Interval gives BMDL and BMDU, which is identical to the
#' confidence interval for the percentile interval approach, SE gives the
#' estimated SE for the "Buckland" method, and modelWeights includes the
#' estimated weights.
#' @author Signe M. Jensen and Jens Riis Baalkilde
#' @references Budtz-Jorgensen, E., Keiding, N., and Grandjean, P. (2001)
#' Benchmark Dose Calculation from Epidemiological Data, \emph{Biometrics}
#' \bold{57}, 698--706.
#' 
#' Crump, K. (2002) Critical Issues in Benchmark Calculations from Continuous
#' Data, \emph{Critical Reviews in Toxicology} \bold{32}, 133--153.
#' 
#' Baalkilde, J. R., Hansen, N. R., and Jensen, S. M. (2025) Stacking Weights
#' and Model Space Selection in Frequentist Model Averaging for Benchmark Dose
#' Estimation, \emph{Environmetrics} \bold{36(2)}, e70002.
#' @keywords model averaging nonlinear bootstrap
#' @export
#' @examples
#' 
#' library(drc)
#' library(drcData)
#' 
#' ## Fitting 4 different two-parameter models to binomial data
#' deguelin.m1 <- drm(r/n~dose, weights=n, data=deguelin, fct=LL.2(), type="binomial")
#' deguelin.m2 <- drm(r/n~dose, weights=n, data=deguelin, fct=W1.2(), type="binomial")
#' deguelin.m3 <- drm(r/n~dose, weights=n, data=deguelin, fct=W2.2(), type="binomial")
#' deguelin.m4 <- drm(r/n~dose, weights=n, data=deguelin, fct=LN.2(), type="binomial")
#' 
#' 
#' ## Model averaged BMD for 5% additional risk with estimated background risk 
#' ## and BMDL based on Buckland et al.
#' bmdMA(list(deguelin.m1,deguelin.m2,deguelin.m3,deguelin.m4), modelWeights="AIC", 0.05, 
#'       backgType = "modelBased", def = "additional",
#'       type = "Buckland")
#' 
#' ## Model averaged BMD for 5% additional risk with estimated background risk
#' ## and BMDL based on an average of the model curves
#' bmdMA(list(deguelin.m1,deguelin.m2,deguelin.m3,deguelin.m4), modelWeights="AIC", 0.05, 
#'       backgType = "modelBased", def = "additional",
#'       type = "curve", bootstrapType = "parametric", bootInterval = "percentile", R=50)
#' 
#' 
#' ## Fitting 4 different two-parameter models to binomial data
#' ryegrass.m1<-drm(rootl~conc, data=ryegrass, fct=LL.4())
#' ryegrass.m2<-drm(rootl~conc, data=ryegrass, fct=W1.4())
#' ryegrass.m3<-drm(rootl~conc, data=ryegrass, fct=W2.4())
#' ryegrass.m4<-drm(rootl~conc, data=ryegrass, fct=LN.4())
#' 
#' ## Model-averaged BMD and bootstrap BMDL for bmr=5% and using the hybrid approach
#' ## to estimate the background risk.  
#' bmdMA(list(ryegrass.m1,ryegrass.m2,ryegrass.m3,ryegrass.m4), modelWeights="AIC", bmr=0.05, 
#'       backgType = "hybridSD", def = "hybridAdd", type = "bootstrap",
#'       bootstrapType = "nonparametric", bootInterval = "percentile", R = 50)
#' 
#' 
#' ## Model-averaged BMD using the Stacking Weights
#' bmdMA(list(ryegrass.m1,ryegrass.m2,ryegrass.m3,ryegrass.m4), modelWeights="Stack", bmr=0.05, 
#'       backgType = "hybridSD", def = "hybridAdd", type = "bootstrap",
#'       bootstrapType = "nonparametric", bootInterval = "percentile", R = 50)
#' 
#' 
bmdMA <- function(modelList, modelWeights, bmr, 
                  backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
                  backg=NA, 
                  def = c("excess", "additional", 
                          "relative", "extra", "added", "hybridExc", "hybridAdd", "point"),
                  respTrans = c("none", "log", "sqrt"),
                  interval = c("delta", "sandwich", "inv", "profile"),
                  type = c("curve","bootstrap","Kang","Buckland"),
                  bootstrapType = c("nonparametric", "parametric"),
                  R=1000,
                  bootInterval = "percentile",
                  level=0.95,
                  stackingSeed = NULL, stackingSplits = 2,
                  display=TRUE, progressInfo = TRUE){
  # assertions
  if(!all(sapply(modelList, function(x) inherits(x, "drc")))){
    stop('modelList must be a list of models of class "drc"')
  }
  
  if(missing(def)) {
    stop(paste("def is missing", sep=""))
  }
  
  if(!modelWeights[1] %in% c("AIC", "BIC", "Stack", "Stacking")){
    if(!(inherits(modelWeights, "numeric") & length(modelWeights) == length(modelList))){
      stop('modelWeights must either be "AIC", "BIC", "Stack", "Stacking" or a numeric vector of same length as modelList')
    }
  }
  
  if(all(type == c("curve","bootstrap","Kang","Buckland")) | sum(type[1] == c("curve","bootstrap", "Bootstrap", "Kang","Buckland")) != 1){
    stop('Specify model averaging type. Options are "curve", "bootstrap", "Kang" and "Buckland"')
  }
  
  if(length(bootstrapType) != 1){
    if(!identical(bootstrapType, c("nonparametric", "parametric"))){
      stop('"bootstrapType" not recognised. Options are: "nonparametric" and "parametric"')
    }
    bootstrapType <- "nonparametric" # default
  }
  if(!bootstrapType %in% c("nonparametric", "parametric")){
    stop('"bootstrapType" not recognised. Options are: "nonparametric" and "parametric"')
  }
  
  interval <- match.arg(interval)
  
  nCurves <- ncol(modelList[[1]]$parmMat)
  bmdList<-lapply(modelList, FUN=function(object){bmd(object, bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans,
                                                      interval = interval, display=FALSE, level=level)})
  
  if(nCurves == 1){
    # Estimate weights
    if(identical(modelWeights,"AIC")){
      modelWeights0<-exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))/2)/
        sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))/2))
    } else if(identical(modelWeights,"BIC")){
      modelWeights0<-exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))/2)/
        sum(exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))/2))
    } else if(identical(modelWeights,"Stack")|identical(modelWeights, "Stacking")){
      # If stackingSeed supplied, save initial seed for later, and set seed for stacking
      if (!is.null(stackingSeed)) {
        sysSeed <- .GlobalEnv$.Random.seed
        set.seed(stackingSeed, kind = "Mersenne-Twister", normal.kind = "Inversion")
      }
      # estimate weights
      modelWeights0 <- getStackingWeights(modelList, stackingSplits)
      
      # If stackingSeed supplied, restore initial seed
      if (!is.null(stackingSeed)) {
        if (!is.null(sysSeed)) {
          .GlobalEnv$.Random.seed <- sysSeed 
        } else {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }
    } else {
      modelWeights0 <- modelWeights
    }
    
    if(identical(modelList[[1]]$type,"continuous")){
      my.fun<-function(x,y){drm(y$call$formula, data = x, fct = y[["fct"]])}
      
      if(identical(type,"Kang")){
        maBMD <- sum(modelWeights0 * sapply(bmdList, function(x){x$Results[,1]}))
        maBMDL <- sum(modelWeights0 * sapply(bmdList, function(x){x$interval[,1]}))
        maBMDU <- sum(modelWeights0 * sapply(bmdList, function(x){x$interval[,2]}))
      }
      
      if(identical(type,"Buckland")){
        if(!(interval %in% c("delta", "sandwich"))){
          stop('Buckland CI can only be estimated for interval type "delta" or "sandwich".')
        }
        estBMD <- sapply(bmdList, function(x){x$Results[,1]})
        seBMD <- sapply(bmdList, function(x){x$SE})
        maBMD <- sum(modelWeights0 * estBMD)
        maBMDse <- sum(modelWeights0 * sqrt(seBMD^2 + (estBMD-maBMD)^2))
        quant <- qnorm(level)
        maBMDL <- maBMD - quant*maBMDse
        maBMDU <- maBMD + quant*maBMDse
      }
      
      if(identical(type,"bootstrap") |identical(type,"Bootstrap") ){
        maBMD <- sum(modelWeights0 * sapply(bmdList, function(x){x$Results[,1]}))
        
        # Bootstrap
        if(identical(bootstrapType,"nonparametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,bootType="nonparametric")
        # } else if(identical(bootstrapType,"semiparametric")){
        #   bootData <- bootDataGen(modelList[[1]],R=R,bootType="semiparametric")
        } else if(identical(bootstrapType,"parametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,bootType="parametric")
        }
        
        bmdMAboot <- function(data){
          bootModelList <- lapply(modelList, function(model) try(
            eval(substitute(drm(formula = formula0, data = data, fct = model$fct, weights = weights0, start = start0,
                                control = drmc(noMessage = TRUE)),
                            list(formula0 = model$call$formula, 
                                 weights0 = model$call$weights,
                                 start0 = coef(model)))),
            silent = TRUE))
          
          modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
          
          bootModelList <- bootModelList[!modelConvergenceError]
          
          if(length(modelWeights) > 1){
            bootModelWeights <- modelWeights[!modelConvergenceError]
          } else {
            bootModelWeights <- modelWeights
          }
          
          bmdEst <- try(bmdMA(bootModelList, bootModelWeights, bmr = bmr, backgType = backgType, 
                              backg = backg, def = def, respTrans = respTrans, interval = "delta", 
                              type = "Kang", stackingSplits = stackingSplits, display = FALSE)$Results[,1], silent = TRUE)
          as.numeric(bmdEst)
        }
        
        bootBmdEst <- numeric(length(bootData))
        
        if(progressInfo){
          cat("Performing bootstrap\n")
          maxIter <- ifelse(bootInterval == "BCa", R + modelList[[1]]$sumList$lenData, R)
          pb <- txtProgressBar(min = 0, max = maxIter, style = 3)
        }
        
        for(i in 1:length(bootData)){
          bootBmdEst[i] <- bmdMAboot(bootData[[i]])
          if(progressInfo) setTxtProgressBar(pb, i)
        }
        if(progressInfo & (bootInterval != "BCa")) close(pb)
        
        boot0<-bootBmdEst[!is.na(bootBmdEst)]
        
        if(length(boot0) == 0){ 
          maBMDL <- NA 
          maBMDU <- NA
        } else {
          if(bootInterval %in% c("percentile","Percentile")){
            maBMDL <- quantile(boot0,p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
            maBMDU <- quantile(boot0,p=c(level), na.rm = TRUE)
          }
          
          if(identical(bootInterval,"BCa")){
            jackData <- list()
            for(i in 1:(dim(modelList[[1]]$data)[1])){
              jackData[[i]] <- modelList[[1]]$data[-i,]
            }
            
            jackBmdEst <- numeric(length(jackData))
            for(i in 1:length(jackData)){
              jackBmdEst[i] <- bmdMAboot(jackData[[i]])
              if(progressInfo) setTxtProgressBar(pb, i + R)
            }
            if(progressInfo) close(pb)
            
            bootjack<-jackBmdEst[!is.na(jackBmdEst)]
            
            maBMDL <- as.numeric(BCa(obs = maBMD, data = modelList[[1]]$data, boot0, bootjack, level = level)[1])
            maBMDU <- "Not available for BCa bootstrap"
          }
        }
      }
      
      if(identical(type,"curve")){
        bmrScaled0<-sum(sapply(bmdList, function(x){x$bmrScaled})*modelWeights0)
        maBMD <- bmdMACurve(modelList,modelWeights0,bmrScaled0)$Results[1]
        
        # Bootstrap
        if(identical(bootstrapType,"nonparametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,bootType="nonparametric")
        # } else if(identical(bootstrapType,"semiparametric")){
        #   bootData <- bootDataGen(modelList[[1]],R=R,bootType="semiparametric")
        } else if(identical(bootstrapType,"parametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,bootType="parametric")
        }
        
        bmdMACurveboot <- function(data){
          bootModelList <- lapply(modelList, function(model){try(
            eval(substitute(
              drm(formula = formula0, data = data, fct = model$fct, weights = weights0, start = start0, type = model$type,
                  control = drmc(noMessage = TRUE)),
              list(formula0 = model$call$formula, 
                   weights0 = model$call$weights,
                   start0 = coef(model)))),
            silent = TRUE)
          })
          
          modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
          
          bootModelList <- bootModelList[!modelConvergenceError]
          bootBmdList <- lapply(bootModelList, 
                                function(object){
                                  try(bmd(object, bmr = bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans,
                                      interval = "delta", display=FALSE), silent = TRUE)})
          
          # Estimate weights
          if(identical(modelWeights,"AIC")){
            bootModelWeights0<-exp(-(sapply(bootModelList,AIC)-min(sapply(bootModelList,AIC))))/
              sum(exp(-(sapply(bootModelList,AIC)-min(sapply(bootModelList,AIC)))))
          } else if(identical(modelWeights,"BIC")){
            bootModelWeights0<-exp(-(sapply(bootModelList,BIC)-min(sapply(bootModelList,BIC))))/
              sum(exp(-(sapply(bootModelList,BIC)-min(sapply(bootModelList,BIC)))))
          } else if(identical(modelWeights,"Stack")|identical(modelWeights, "Stacking")){
            # estimate weights
            bootModelWeights0 <- getStackingWeights(bootModelList, stackingSplits)
          } else {
            bootModelWeights0 <- modelWeights[!modelConvergenceError]
          }
          
          # bootBmrScaled0 <- colSums(bootModelWeights0 * t(sapply(bootBmdList, function(x){x$bmrScaled})))
          bootBmrScaled0 <- sum(sapply(bootBmdList, function(x){as.numeric(try(x$bmrScaled, TRUE))})*bootModelWeights0)
          
          # bootBmdEst <- try(bmdMACurve(bootModelList,bootModelWeights0,bootBmrScaled0)$Results[,1], silent = TRUE)
          bootBmdEst <- try(bmdMACurve(bootModelList,bootModelWeights0,bootBmrScaled0)$Results[1], silent = TRUE)
          as.numeric(bootBmdEst)
        }
        
        bootBmdEst <- numeric(length(bootData))
        
        if(progressInfo){
          cat("Performing bootstrap\n")
          maxIter <- ifelse(bootInterval == "BCa", R + modelList[[1]]$sumList$lenData, R)
          pb <- txtProgressBar(min = 0, max = maxIter, style = 3)
        }
        
        for(i in 1:length(bootData)){
          bootBmdEst[i] <- bmdMACurveboot(bootData[[i]])
          if(progressInfo) setTxtProgressBar(pb, i)
        }
        if(progressInfo & (bootInterval != "BCa")) close(pb)
        
        boot0 <- bootBmdEst[!is.na(bootBmdEst)]
        
        if(length(boot0) == 0){ 
          maBMDL <- NA 
          maBMDU <- NA
        } else {
          if(bootInterval %in% c("percentile","Percentile")){
            maBMDL <- quantile(boot0,p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
            maBMDU <- quantile(boot0,p=c(level), na.rm = TRUE)
          }
          
          if(identical(bootInterval,"BCa")){
            jackData <- list()
            for(i in 1:(dim(modelList[[1]]$data)[1])){
              jackData[[i]] <- modelList[[1]]$data[-i,]
            }
            
            jackBmdEst <- numeric(length(jackData))
            
            for(i in 1:length(jackData)){
              jackBmdEst[i] <- bmdMACurveboot(jackData[[i]])
              if(progressInfo) setTxtProgressBar(pb, i + R)
            }
            if(progressInfo) close(pb)
            
            bootjack <- jackBmdEst[!is.na(jackBmdEst)]
            
            maBMDL <- as.numeric(BCa(obs = maBMD, data = modelList[[1]]$data, boot0, bootjack, level = level)[1])
            maBMDU <- "Not available for BCa bootstrap"
          }
        }
      }
    }
    
    if(identical(modelList[[1]]$type,"binomial")){
      my.fun<-function(x,y){drm(y$call$formula, data = x, type="binomial", fct = y[["fct"]])}
      my.fun2<-function(x,y){drm(number~dose, data = x, type="binomial", fct = y[["fct"]])}
      
      if(identical(type,"Kang")){
        maBMD <- sum(modelWeights0 * sapply(bmdList, function(x){x$Results[,1]}))
        maBMDL <- sum(modelWeights0 * sapply(bmdList, function(x){x$interval[,1]}))
        maBMDU <- sum(modelWeights0 * sapply(bmdList, function(x){x$interval[,2]}))
      }
      
      if(identical(type,"Buckland")){
        if(!(interval %in% c("delta", "sandwich"))){
          stop('Buckland CI can only be estimated for interval type "delta" or "sandwich".')
        }
        estBMD <- sapply(bmdList, function(x){x$Results[,1]})
        seBMD <- sapply(bmdList, function(x){x$SE})
        maBMD <- sum(modelWeights0 * estBMD)
        maBMDse <- sum(modelWeights0 * sqrt(seBMD^2 + (estBMD-maBMD)^2))
        quant <- qnorm(level)
        maBMDL <- maBMD - quant*maBMDse
        maBMDU <- maBMD + quant*maBMDse
      }
      
      if(identical(type,"bootstrap") | identical(type,"Bootstrap") ){
        maBMD <- sum(modelWeights0 * sapply(bmdList, function(x){x$Results[,1]}))
        
        if(identical(bootstrapType,"nonparametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,bootType="nonparametric",aggregated = FALSE)
        # } else if(identical(bootstrapType,"semiparametric")){
        #   bootData <- bootDataGen(modelList[[1]],R=R,bootType="semiparametric",aggregated = FALSE)
        } else if(identical(bootstrapType,"parametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,bootType="parametric",aggregated = FALSE)
        }
        
        bmdMAboot <- function(data){
          bootModelList <- lapply(modelList, function(model) try(
            eval(substitute(drm(formula = formula0, data = data, fct = model$fct, weights = weights0, start = start0, type = "binomial",
                                control = drmc(noMessage = TRUE)),
                            list(formula0 = model$call$formula, 
                                 weights0 = model$call$weights,
                                 start0 = coef(model)))),
            silent = TRUE))
          
          modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
          
          bootModelList <- bootModelList[!modelConvergenceError]
          
          if(length(modelWeights) > 1){
            bootModelWeights <- modelWeights[!modelConvergenceError]
          } else {
            bootModelWeights <- modelWeights
          }
          
          bmdEst <- try(bmdMA(bootModelList, bootModelWeights, bmr = bmr, backgType = backgType, 
                              backg = backg, def = def, respTrans = respTrans, interval = "delta", 
                              type = "Kang", stackingSplits = stackingSplits, display = FALSE)$Results[,1], silent = TRUE)
          as.numeric(bmdEst)
        }
        
        bootBmdEst <- numeric(length(bootData))
        
        if(progressInfo){
          cat("Performing bootstrap\n")
          data.str <- modelList[[1]]$data
          maxIter <- ifelse(bootInterval == "BCa", R + sum(data.str[["weights"]]), R)
          pb <- txtProgressBar(min = 0, max = maxIter, style = 3)
        }
        
        for(i in 1:length(bootData)){
          bootBmdEst[i] <- bmdMAboot(bootData[[i]])
          if(progressInfo) setTxtProgressBar(pb, i)
        }
        if(progressInfo & (bootInterval != "BCa")) close(pb)
        
        boot0<-bootBmdEst[!is.na(bootBmdEst)]
        
        if(length(boot0) == 0){ 
          maBMDL <- NA 
          maBMDU <- NA
        } else {
          if(bootInterval %in% c("percentile","Percentile")){
            maBMDL <- quantile(boot0,p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
            maBMDU <- quantile(boot0,p=c(level), na.rm = TRUE)
          }
          
          if(identical(bootInterval,"BCa")){
            data.str <- modelList[[1]]$data
            data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
            data.e<-expandBinomial(data.str, 
                                   number = "number",
                                   total = "weights",
                                   dose = as.character(modelList[[1]]$call$formula[[3]]))
            df <- data.frame(data.e[,as.character(modelList[[1]]$call$formula[[3]])],
                             data.e[,"number"],
                             data.e[,"weights"])
            colnames(df) <- c(as.character(modelList[[1]]$call$formula[[3]]),
                              as.character(modelList[[1]]$call$formula[[2]])[[2]],
                              as.character(modelList[[1]]$call$formula[[2]])[[3]])
            jackData <- list()
            for(i in 1:(dim(df)[1])){
              jackData[[i]] <- df[-i,]
            }
            
            jackBmdEst <- numeric(length(jackData))
            for(i in 1:length(jackData)){
              jackBmdEst[i] <- bmdMAboot(jackData[[i]])
              if(progressInfo) setTxtProgressBar(pb, i + R)
            }
            if(progressInfo) close(pb)
            
            bootjack<-jackBmdEst[!is.na(jackBmdEst)]
            
            maBMDL <- as.numeric(BCa(obs = maBMD, data = df, boot0, bootjack, level = level)[1]) # data = modelList[[1]]$data
            maBMDU <- "Not available for BCa bootstrap"
          }
        }
      }
      
      if(identical(type,"curve")){
        bmrScaled0<-sum(sapply(bmdList, function(x){x$bmrScaled})*modelWeights0)
        maBMD <- bmdMACurve(modelList,modelWeights0,bmrScaled0)$Results[1]
        
        # Bootstrap
        if(identical(bootstrapType,"nonparametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,bootType="nonparametric",aggregated = FALSE)
        # } else if(identical(bootstrapType,"semiparametric")){
        #   bootData <- bootDataGen(modelList[[1]],R=R,bootType="semiparametric",aggregated=FALSE)
        } else if(identical(bootstrapType,"parametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,bootType="parametric",aggregated=FALSE)
        }
        
        bmdMACurveboot <- function(data){
          bootModelList <- lapply(modelList, function(model){try(
            eval(substitute(
              drm(formula = formula0, data = data, fct = model$fct, weights = weights0, start = start0, type = model$type,
                  control = drmc(noMessage = TRUE)),
              list(formula0 = model$call$formula, 
                   weights0 = model$call$weights,
                   start0 = coef(model)))),
            silent = TRUE)
          })
          
          modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
          
          bootModelList <- bootModelList[!modelConvergenceError]
          bootBmdList <- lapply(bootModelList, 
                                function(object){
                                  try(bmd(object, bmr = bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans,
                                      interval = "delta", display=FALSE), silent = TRUE)})
          
          # Estimate weights
          if(identical(modelWeights,"AIC")){
            bootModelWeights0<-exp(-(sapply(bootModelList,AIC)-min(sapply(bootModelList,AIC))))/
              sum(exp(-(sapply(bootModelList,AIC)-min(sapply(bootModelList,AIC)))))
          } else if(identical(modelWeights,"BIC")){
            bootModelWeights0<-exp(-(sapply(bootModelList,BIC)-min(sapply(bootModelList,BIC))))/
              sum(exp(-(sapply(bootModelList,BIC)-min(sapply(bootModelList,BIC)))))
          } else if(identical(modelWeights,"Stack")|identical(modelWeights, "Stacking")){
            # estimate weights
            bootModelWeights0 <- getStackingWeights(bootModelList, stackingSplits)
          } else {
            bootModelWeights0 <- modelWeights[!modelConvergenceError]
          }
          
          
          bootBmrScaled0 <- sum(sapply(bootBmdList, function(x){as.numeric(try(x$bmrScaled, TRUE))})*bootModelWeights0)
          
          bootBmdEst <- try(bmdMACurve(bootModelList,bootModelWeights0,bootBmrScaled0)$Results[1], silent = TRUE)
          as.numeric(bootBmdEst)
        }
        
        bootBmdEst <- numeric(length(bootData))
        
        if(progressInfo){
          cat("Performing bootstrap\n")
          data.str <- modelList[[1]]$data
          maxIter <- ifelse(bootInterval == "BCa", R + sum(data.str[["weights"]]), R)
          pb <- txtProgressBar(min = 0, max = maxIter, style = 3)
        }
        
        for(i in 1:length(bootData)){
          bootBmdEst[i] <- bmdMACurveboot(bootData[[i]])
          if(progressInfo) setTxtProgressBar(pb, i)
        }
        if(progressInfo & (bootInterval != "BCa")) close(pb)
        
        boot0 <- bootBmdEst[!is.na(bootBmdEst)]
        
        if(length(boot0) == 0){ 
          maBMDL <- NA 
          maBMDU <- NA
        } else {
          if(bootInterval %in% c("percentile","Percentile")){
            maBMDL <- quantile(boot0,p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
            maBMDU <- quantile(boot0,p=c(level), na.rm = TRUE)
          }
          
          if(identical(bootInterval,"BCa")){
            data.str <- modelList[[1]]$data
            data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
            data.e<-expandBinomial(data.str, 
                                   number = "number",
                                   total = "weights",
                                   dose = as.character(modelList[[1]]$call$formula[[3]]))
            df <- data.frame(data.e[,as.character(modelList[[1]]$call$formula[[3]])],
                             data.e[,"number"],
                             data.e[,"weights"])
            colnames(df) <- c(as.character(modelList[[1]]$call$formula[[3]]),
                              as.character(modelList[[1]]$call$formula[[2]])[[2]],
                              as.character(modelList[[1]]$call$formula[[2]])[[3]])
            jackData <- list()
            for(i in 1:nrow(df)){
              jackData[[i]] <- df[-i,]
            }
            
            jackBmdEst <- numeric(length(jackData))
            
            for(i in 1:length(jackData)){
              jackBmdEst[i] <- bmdMACurveboot(jackData[[i]])
              if(progressInfo) setTxtProgressBar(pb, i + R)
            }
            if(progressInfo) close(pb)
            
            bootjack <- jackBmdEst[!is.na(jackBmdEst)]
            
            maBMDL <- as.numeric(BCa(obs = maBMD, data = df, boot0, bootjack, level = level)[1]) # data = modelList[[1]]$data
            maBMDU <- "Not available for BCa bootstrap"
          }
        }
      }
    }
  } 
  
  if (nCurves > 1){
    if(is.null(modelList[[1]]$objList)){
      # Estimate weights
      if(identical(modelWeights,"AIC")){
        modelWeights0<-exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))/2)/
          sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))/2))
      } else if(identical(modelWeights,"BIC")){
        modelWeights0<-exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))/2)/
          sum(exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))/2))
      } else if(identical(modelWeights,"Stack")|identical(modelWeights, "Stacking")){
        # If stackingSeed supplied, save initial seed for later, and set seed for stacking
        if (!is.null(stackingSeed)) {
          sysSeed <- .GlobalEnv$.Random.seed
          set.seed(stackingSeed, kind = "Mersenne-Twister", normal.kind = "Inversion")
        }
        # estimate weights
        modelWeights0 <- getStackingWeights(modelList, stackingSplits)
        
        # If stackingSeed supplied, restore initial seed
        if (!is.null(stackingSeed)) {
          if (!is.null(sysSeed)) {
            .GlobalEnv$.Random.seed <- sysSeed 
          } else {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        }
      } else {
        modelWeights0 <- modelWeights
      }
      
      if(identical(modelList[[1]]$type,"continuous")){
        if(identical(type,"Kang")){
          maBMD <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$Results[,1])))
          maBMDL <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$interval[,1])))
          maBMDU <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$interval[,2])))
        }
        
        if(identical(type,"Buckland")){
          if(!(interval %in% c("delta", "sandwich"))){
            stop('Buckland CI can only be estimated for interval type "delta" or "sandwich".')
          }
          estBMD <- sapply(bmdList, function(x) x$Results[,1]) #sapply(bmdList, function(x){x$Results[,1]})
          seBMD <- sapply(bmdList, function(x){x$SE})
          maBMD <- colSums(modelWeights0 * t(estBMD))
          maBMDse <- colSums(modelWeights0 * t(sqrt(seBMD^2 + (estBMD-maBMD)^2)))
          quant <- qnorm(level)
          maBMDL <- maBMD - quant*maBMDse
          maBMDU <- maBMD + quant*maBMDse
        }
        
        if(identical(type,"bootstrap") |identical(type,"Bootstrap") ){
          # Estimate BMD
          maBMD <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$Results[,1])))
          
          # Bootstrap
          if(identical(bootstrapType,"nonparametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,bootType="nonparametric")
          # } else if(identical(bootstrapType,"semiparametric")){
          #   bootData <- bootDataGen(modelList[[1]],R=R,bootType="semiparametric")
          } else if(identical(bootstrapType,"parametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,bootType="parametric")
          }
          
          bmdMAboot <- function(data){
            # data[[as.character(modelList[[1]]$call$curveid)]] <- data[[paste0("orig.", as.character(modelList[[1]]$call$curveid))]]
            if(is.null(modelList[[1]]$call$pmodels)){
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula = formula0, data = data, fct = model$fct, weights = weights0, start = start0, 
                      curveid = curveid0, type = model$type, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula, 
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       start0 = coef(model)))),
                silent = TRUE)
              })
            } else {
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula0, weights = weights0, curveid = curveid0, pmodels = pmodels0, start = start0,
                      data = data, type = model$type, fct = model$fct, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula,
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       pmodels0 = model$call$pmodels,
                       start0 = coef(model))
                )),
                silent = TRUE)
              })
            }
            modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
            
            bootModelList <- bootModelList[!modelConvergenceError]
            
            if(length(modelWeights) > 1){
              bootModelWeights <- modelWeights[!modelConvergenceError]
            } else {
              bootModelWeights <- modelWeights
            }
            
            bmdEst <- try(bmdMA(bootModelList, bootModelWeights, bmr = bmr, backgType = backgType, 
                                backg = backg, def = def, respTrans = respTrans, interval = "delta", 
                                type = "Kang", stackingSplits = stackingSplits, display = FALSE)$Results[colnames(modelList[[1]]$parmMat),1], silent = TRUE)
            as.numeric(bmdEst)
          }
          
          bootBmdEst <- matrix(NA, nrow = length(bootData), ncol = ncol(modelList[[1]]$parmMat))
          
          if(progressInfo){
            cat("Performing bootstrap\n")
            maxIter <- ifelse(bootInterval == "BCa", R + modelList[[1]]$sumList$lenData, R)
            pb <- txtProgressBar(min = 0, max = maxIter, style = 3)
          }
          
          for(i in 1:length(bootData)){
            bootBmdEst[i,] <- bmdMAboot(bootData[[i]])
            if(progressInfo) setTxtProgressBar(pb, i)
          }
          if(progressInfo & (bootInterval != "BCa")) close(pb)
          
          bootBmdError <- apply(bootBmdEst, 1, function(x) any(is.na(x)))
          boot0 <- bootBmdEst[!bootBmdError,]
          
          if(length(boot0) == 0){ 
            maBMDL <- NA 
            maBMDU <- NA
          } else {
            if(bootInterval %in% c("percentile","Percentile")){
              maBMDL <- apply(boot0, 2, quantile, p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
              maBMDU <- apply(boot0, 2, quantile, p=c(level), na.rm = TRUE)
            }
            
            if(identical(bootInterval,"BCa")){
              jackData <- list()
              for(i in 1:(dim(modelList[[1]]$data)[1])){
                jackData[[i]] <- modelList[[1]]$origData[-i,]
              }
              
              jackBmdEst <- matrix(NA, nrow = length(jackData), ncol = nCurves)
              for(i in 1:length(jackData)){
                jackBmdEst[i,] <- bmdMAboot(jackData[[i]])
                if(progressInfo) setTxtProgressBar(pb, i + R)
              }
              if(progressInfo) close(pb)
              
              jackBmdError <- apply(jackBmdEst, 1, function(x) any(is.na(x)))
              bootjack <- jackBmdEst[!jackBmdError,]
              
              maBMDL <- sapply(1:nCurves, function(i) as.numeric(BCa(obs = maBMD[i], data = modelList[[1]]$data, boot0[,i], bootjack[,i], level = level)[1]))
              maBMDU <- rep("Not available for BCa bootstrap", nCurves)
            }
          }
        }
        
        if(identical(type,"curve")){
          bmrScaled0 <- colSums(modelWeights0 * t(sapply(bmdList, function(x){x$bmrScaled})))
          
          maBMD <- bmdMACurve(modelList,modelWeights0,bmrScaled0)$Results[,1]
          
          # Bootstrap
          if(identical(bootstrapType,"nonparametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,bootType="nonparametric")
          # } else if(identical(bootstrapType,"semiparametric")){
          #   bootData <- bootDataGen(modelList[[1]],R=R,bootType="semiparametric")
          } else if(identical(bootstrapType,"parametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,bootType="parametric")
          }
          
          bmdMACurveboot <- function(data){
            # data[[as.character(modelList[[1]]$call$curveid)]] <- data[[paste0("orig.", as.character(modelList[[1]]$call$curveid))]]
            if(is.null(modelList[[1]]$call$pmodels)){
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula = formula0, data = data, fct = model$fct, weights = weights0, start = start0, 
                      curveid = curveid0, type = model$type, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula, 
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       start0 = coef(model)))),
                silent = TRUE)
              })
            } else {
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula0, weights = weights0, curveid = curveid0, pmodels = pmodels0, start = start0,
                      data = data, type = model$type, fct = model$fct, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula,
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       pmodels0 = model$call$pmodels,
                       start0 = coef(model))
                )),
                silent = TRUE)
              })
            }
            modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
            
            bootModelList <- bootModelList[!modelConvergenceError]
            bootBmdList <- lapply(bootModelList, 
                                  function(object){
                                    try(bmd(object, bmr = bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans,
                                        interval = "delta", display=FALSE), silent = TRUE)})
            
            # Estimate weights
            if(identical(modelWeights,"AIC")){
              bootModelWeights0<-exp(-(sapply(bootModelList,AIC)-min(sapply(bootModelList,AIC))))/
                sum(exp(-(sapply(bootModelList,AIC)-min(sapply(bootModelList,AIC)))))
            } else if(identical(modelWeights,"BIC")){
              bootModelWeights0<-exp(-(sapply(bootModelList,BIC)-min(sapply(bootModelList,BIC))))/
                sum(exp(-(sapply(bootModelList,BIC)-min(sapply(bootModelList,BIC)))))
            } else if(identical(modelWeights,"Stack")|identical(modelWeights, "Stacking")){
              # estimate weights
              bootModelWeights0 <- getStackingWeights(bootModelList, stackingSplits)
            } else {
              bootModelWeights0 <- modelWeights[!modelConvergenceError]
            }
            
            bootBmrScaled0 <- try(colSums(bootModelWeights0 * t(sapply(bootBmdList, function(x){x$bmrScaled[colnames(modelList[[1]]$parmMat),1]}))), silent = TRUE)
            
            bootBmdEst <- try(bmdMACurve(bootModelList,bootModelWeights0,bootBmrScaled0)$Results[colnames(modelList[[1]]$parmMat),1], silent = TRUE)
            as.numeric(bootBmdEst)
          }
          
          bootBmdEst <- matrix(NA, nrow = length(bootData), ncol = ncol(modelList[[1]]$parmMat))
          
          if(progressInfo){
            cat("Performing bootstrap\n")
            maxIter <- ifelse(bootInterval == "BCa", R + modelList[[1]]$sumList$lenData, R)
            pb <- txtProgressBar(min = 0, max = maxIter, style = 3)
          }
          
          for(i in 1:length(bootData)){
            bootBmdEst[i,] <- bmdMACurveboot(bootData[[i]])
            if(progressInfo) setTxtProgressBar(pb, i)
          }
          if(progressInfo & (bootInterval != "BCa")) close(pb)
          
          bootBmdError <- apply(bootBmdEst, 1, function(x) any(is.na(x)))
          boot0 <- bootBmdEst[!bootBmdError,]
          
          if(length(boot0) == 0){ 
            maBMDL <- NA 
            maBMDU <- NA
          } else {
            if(bootInterval %in% c("percentile","Percentile")){
              maBMDL <- apply(boot0, 2, quantile, p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
              maBMDU <- apply(boot0, 2, quantile, p=c(level), na.rm = TRUE)
            }
            
            if(identical(bootInterval,"BCa")){
              jackData <- list()
              for(i in 1:(dim(modelList[[1]]$data)[1])){
                jackData[[i]] <- modelList[[1]]$origData[-i,]
              }
              
              jackBmdEst <- matrix(NA, nrow = length(jackData), ncol = nCurves)
              
              for(i in 1:length(jackData)){
                jackBmdEst[i,] <- bmdMACurveboot(jackData[[i]])
                if(progressInfo) setTxtProgressBar(pb, i + R)
              }
              if(progressInfo) close(pb)
              
              
              jackBmdError <- apply(jackBmdEst, 1, function(x) any(is.na(x)))
              bootjack <- jackBmdEst[!jackBmdError,]
              
              maBMDL <- sapply(1:nCurves, function(i) as.numeric(BCa(obs = maBMD[i], data = modelList[[1]]$data, boot0[,i], bootjack[,i], level = level)[1]))
              maBMDU <- rep("Not available for BCa bootstrap", nCurves)
            }
          }
        }
      }
      
      if(identical(modelList[[1]]$type,"binomial")){
        if(identical(type,"Kang")){
          maBMD <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$Results[,1])))
          maBMDL <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$interval[,1])))
          maBMDU <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$interval[,2])))
        }
        
        if(identical(type,"Buckland")){
          if(!(interval %in% c("delta", "sandwich"))){
            stop('Buckland CI can only be estimated for interval type "delta" or "sandwich".')
          }
          estBMD <- sapply(bmdList, function(x) x$Results[,1]) #sapply(bmdList, function(x){x$Results[,1]})
          seBMD <- sapply(bmdList, function(x){x$SE})
          maBMD <- colSums(modelWeights0 * t(estBMD))
          maBMDse <- colSums(modelWeights0 * t(sqrt(seBMD^2 + (estBMD-maBMD)^2)))
          quant <- qnorm(level)
          maBMDL <- maBMD - quant*maBMDse
          maBMDU <- maBMD + quant*maBMDse
        }
        
        if(identical(type,"bootstrap") |identical(type,"Bootstrap") ){
          # Estimate BMD
          maBMD <- colSums(modelWeights0 * t(sapply(bmdList, function(x) x$Results[,1])))
          
          # Bootstrap
          if(identical(bootstrapType,"nonparametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,bootType="nonparametric",aggregated = FALSE)
          # } else if(identical(bootstrapType,"semiparametric")){
          #   bootData <- bootDataGen(modelList[[1]],R=R,bootType="semiparametric",aggregated=FALSE)
          } else if(identical(bootstrapType,"parametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,bootType="parametric",aggregated=FALSE)
          }
          
          bmdMAboot <- function(data){
            if(is.null(modelList[[1]]$call$pmodels)){
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula = formula0, data = data, fct = model$fct, weights = weights0, start = start0, 
                      curveid = curveid0, type = model$type, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula, 
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       start0 = coef(model)))),
                silent = TRUE)
              })
            } else {
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula0, weights = weights0, curveid = curveid0, pmodels = pmodels0, start = start0,
                      data = data, type = model$type, fct = model$fct, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula,
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       pmodels0 = model$call$pmodels,
                       start0 = coef(model))
                )),
                silent = TRUE)
              })
            }
            modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
            
            bootModelList <- bootModelList[!modelConvergenceError]
            
            if(length(modelWeights) > 1){
              bootModelWeights <- modelWeights[!modelConvergenceError]
            } else {
              bootModelWeights <- modelWeights
            }
            
            bmdEst <- try(bmdMA(bootModelList, bootModelWeights, bmr = bmr, backgType = backgType, 
                                backg = backg, def = def, respTrans = respTrans, interval = "delta", 
                                type = "Kang", stackingSplits = stackingSplits, display = FALSE)$Results[colnames(modelList[[1]]$parmMat),1], silent = TRUE)
            as.numeric(bmdEst)
          }
          
          bootBmdEst <- matrix(NA, nrow = length(bootData), ncol = ncol(modelList[[1]]$parmMat))
          
          if(progressInfo){
            cat("Performing bootstrap\n")
            data.str <- modelList[[1]]$data
            maxIter <- ifelse(bootInterval == "BCa", R + sum(data.str[["weights"]]), R)
            pb <- txtProgressBar(min = 0, max = maxIter, style = 3)
          }
          
          for(i in 1:length(bootData)){
            bootBmdEst[i,] <- bmdMAboot(bootData[[i]])
            if(progressInfo) setTxtProgressBar(pb, i)
          }
          if(progressInfo & (bootInterval != "BCa")) close(pb)
          
          bootBmdError <- apply(bootBmdEst, 1, function(x) any(is.na(x)))
          boot0 <- bootBmdEst[!bootBmdError,]
          
          if(length(boot0) == 0){ 
            maBMDL <- NA 
            maBMDU <- NA
          } else {
            if(bootInterval %in% c("percentile","Percentile")){
              maBMDL <- apply(boot0, 2, quantile, p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
              maBMDU <- apply(boot0, 2, quantile, p=c(level), na.rm = TRUE)
            }
            
            if(identical(bootInterval,"BCa")){
              data.str <- modelList[[1]]$data
              data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
              data.e<-expandBinomial(data.str, 
                                     number = "number",
                                     total = "weights",
                                     dose = as.character(modelList[[1]]$call$formula[[3]]),
                                     curveid = as.character(modelList[[1]]$call$curveid))
              df <- data.frame(data.e[,as.character(modelList[[1]]$call$formula[[3]])],
                               data.e[,"number"],
                               data.e[,"weights"],
                               data.e[,as.character(modelList[[1]]$call$curveid)])
              colnames(df) <- c(as.character(modelList[[1]]$call$formula[[3]]),
                                as.character(modelList[[1]]$call$formula[[2]])[[2]],
                                as.character(modelList[[1]]$call$formula[[2]])[[3]],
                                as.character(modelList[[1]]$call$curveid))
              
              jackData <- list()
              for(i in 1:nrow(df)){
                jackData[[i]] <- df[-i,]
              }
              
              jackBmdEst <- matrix(NA, nrow = length(jackData), ncol = nCurves)
              for(i in 1:length(jackData)){
                jackBmdEst[i,] <- bmdMAboot(jackData[[i]])
                if(progressInfo) setTxtProgressBar(pb, i + R)
              }
              if(progressInfo) close(pb)
              
              jackBmdError <- apply(jackBmdEst, 1, function(x) any(is.na(x)))
              bootjack <- jackBmdEst[!jackBmdError,]
              
              maBMDL <- sapply(1:nCurves, function(i) as.numeric(BCa(obs = maBMD[i], data = df, boot0[,i], bootjack[,i], level = level)[1])) #  data = modelList[[1]]$data
              maBMDU <- rep("Not available for BCa bootstrap", nCurves)
            }
          }
        }
        
        if(identical(type,"curve")){
          bmrScaled0 <- colSums(modelWeights0 * t(sapply(bmdList, function(x){x$bmrScaled})))
          
          maBMD <- bmdMACurve(modelList,modelWeights0,bmrScaled0)$Results[,1]
          
          # Bootstrap
          if(identical(bootstrapType,"nonparametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,bootType="nonparametric",aggregated = FALSE)
          # } else if(identical(bootstrapType,"semiparametric")){
          #   bootData <- bootDataGen(modelList[[1]],R=R,bootType="semiparametric",aggregated=FALSE)
          } else if(identical(bootstrapType,"parametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,bootType="parametric",aggregated=FALSE)
          }
          
          bmdMACurveboot <- function(data){
            if(is.null(modelList[[1]]$call$pmodels)){
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula = formula0, data = data, fct = model$fct, weights = weights0, start = start0,
                      curveid = curveid0, type = model$type, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula, 
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       start0 = coef(model)))),
                silent = TRUE)
              })
            } else {
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula0, weights = weights0, curveid = curveid0, pmodels = pmodels0, start = start0,
                      data = data, type = model$type, fct = model$fct, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula,
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       pmodels0 = model$call$pmodels,
                       start0 = coef(model))
                )),
                silent = TRUE)
              })
            }
            modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
            
            bootModelList <- bootModelList[!modelConvergenceError]
            bootBmdList <- lapply(bootModelList, 
                                  function(object){
                                    try(bmd(object, bmr = bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans,
                                        interval = "delta", display=FALSE), silent = TRUE)})
            
            # Estimate weights
            if(identical(modelWeights,"AIC")){
              bootModelWeights0<-exp(-(sapply(bootModelList,AIC)-min(sapply(bootModelList,AIC))))/
                sum(exp(-(sapply(bootModelList,AIC)-min(sapply(bootModelList,AIC)))))
            } else if(identical(modelWeights,"BIC")){
              bootModelWeights0<-exp(-(sapply(bootModelList,BIC)-min(sapply(bootModelList,BIC))))/
                sum(exp(-(sapply(bootModelList,BIC)-min(sapply(bootModelList,BIC)))))
            } else if(identical(modelWeights,"Stack")|identical(modelWeights, "Stacking")){
              # estimate weights
              bootModelWeights0 <- getStackingWeights(bootModelList, stackingSplits)
            } else {
              bootModelWeights0 <- modelWeights[!modelConvergenceError]
            }
            
            bootBmrScaled0 <- colSums(bootModelWeights0 * t(sapply(bootBmdList, function(x){x$bmrScaled[colnames(modelList[[1]]$parmMat),1]})))
            
            bootBmdEst <- try(bmdMACurve(bootModelList,bootModelWeights0,bootBmrScaled0)$Results[colnames(modelList[[1]]$parmMat),1], silent = TRUE)
            as.numeric(bootBmdEst)
          }
          
          bootBmdEst <- matrix(NA, nrow = length(bootData), ncol = ncol(modelList[[1]]$parmMat))
          
          if(progressInfo){
            cat("Performing bootstrap\n")
            data.str <- modelList[[1]]$data
            maxIter <- ifelse(bootInterval == "BCa", R + sum(data.str[["weights"]]), R)
            pb <- txtProgressBar(min = 0, max = maxIter, style = 3)
          }
          
          for(i in 1:length(bootData)){
            bootBmdEst[i,] <- bmdMACurveboot(bootData[[i]])
            if(progressInfo) setTxtProgressBar(pb, i)
          }
          if(progressInfo & (bootInterval != "BCa")) close(pb)
          
          bootBmdError <- apply(bootBmdEst, 1, function(x) any(is.na(x)))
          boot0 <- bootBmdEst[!bootBmdError,]
          
          if(length(boot0) == 0){ 
            maBMDL <- NA 
            maBMDU <- NA
          } else {
            if(bootInterval %in% c("percentile","Percentile")){
              maBMDL <- apply(boot0, 2, quantile, p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
              maBMDU <- apply(boot0, 2, quantile, p=c(level), na.rm = TRUE)
            }
            
            if(identical(bootInterval,"BCa")){
              data.str <- modelList[[1]]$data
              data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
              data.e<-expandBinomial(data.str, 
                                     number = "number",
                                     total = "weights",
                                     dose = as.character(modelList[[1]]$call$formula[[3]]),
                                     curveid = as.character(modelList[[1]]$call$curveid))
              df <- data.frame(data.e[,as.character(modelList[[1]]$call$formula[[3]])],
                               data.e[,"number"],
                               data.e[,"weights"],
                               data.e[,as.character(modelList[[1]]$call$curveid)])
              colnames(df) <- c(as.character(modelList[[1]]$call$formula[[3]]),
                                as.character(modelList[[1]]$call$formula[[2]])[[2]],
                                as.character(modelList[[1]]$call$formula[[2]])[[3]],
                                as.character(modelList[[1]]$call$curveid))
              
              jackData <- list()
              for(i in 1:nrow(df)){
                jackData[[i]] <- df[-i,]
              }
              
              jackBmdEst <- matrix(NA, nrow = length(jackData), ncol = nCurves)
              
              for(i in 1:length(jackData)){
                jackBmdEst[i,] <- bmdMACurveboot(jackData[[i]])
                if(progressInfo) setTxtProgressBar(pb, i + R)
              }
              if(progressInfo) close(pb)
              
              
              jackBmdError <- apply(jackBmdEst, 1, function(x) any(is.na(x)))
              bootjack <- jackBmdEst[!jackBmdError,]
              
              maBMDL <- sapply(1:nCurves, function(i) as.numeric(BCa(obs = maBMD[i], data = df, boot0[,i], bootjack[,i], level = level)[1])) # data = modelList[[1]]$data
              maBMDU <- rep("Not available for BCa bootstrap", nCurves)
            }
          }
        }
      }
    } else {
      # CURVES FITTED INDEPENDENTLY
      modelListList <- lapply(1:length(modelList[[1]]$objList), function(i) lapply(modelList, function(object) object$objList[[i]]))
      
      bmdMACall <- function(modelList){
        bmdMA(modelList, modelWeights, bmr, backgType, backg, def, respTrans, 
              interval, type, bootstrapType, R, bootInterval, level, stackingSeed, stackingSplits, display = FALSE)
      }
      
      bmdMAList <- lapply(modelListList, bmdMACall)
      
      maBMD <- sapply(bmdMAList, function(x) x$Results[,1])
      maBMDL <- sapply(bmdMAList, function(x) x$interval[,1])
      maBMDU <- sapply(bmdMAList, function(x) x$interval[,2])
      maBMDse <- sapply(bmdMAList, function(x) x$SE[,1])
      
      modelWeights0 <- do.call(rbind, lapply(bmdMAList, function(x) x$modelWeights))
      rownames(modelWeights0) <- names(modelList[[1]]$objList)
    }
  }
  
  resMat<-matrix(c(maBMD,maBMDL),nCurves,2, byrow = FALSE)
  colnames(resMat) <- c("BMD_MA", "BMDL_MA")
  
  intMat<-matrix(c(maBMDL,maBMDU),nCurves,2, byrow = FALSE)
  colnames(intMat) <- c("BMDL_MA", "BMDU_MA")
  
  if(!identical(type, "Buckland")){
    maBMDse <- NA
  }
  seMat <- matrix(maBMDse, nrow = nCurves)
  colnames(seMat) <- "SE"
  
  if(nCurves == 1){
    rownames(resMat) <- ""
    rownames(intMat) <- ""
    rownames(seMat) <- ""
  } else {
    rownames(resMat) <- colnames(modelList[[1]]$parmMat)
    rownames(intMat) <- colnames(modelList[[1]]$parmMat)
    rownames(seMat) <- colnames(modelList[[1]]$parmMat)
  }
  
  used.Boot<-ifelse(identical(type,"bootstrap")|identical(type,"Bootstrap")|identical(type,"curve"),
                    floor(length(boot0)/nCurves),NA)
  resBMD<-list(Results = resMat,
               Boot.samples.used = used.Boot,
               interval = intMat,
               SE = seMat,
               modelWeights = modelWeights0)
  if(display){
    print(resMat)
  }
  
  class(resBMD) <- "bmd"
  invisible(resBMD)
}




