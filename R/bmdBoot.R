#' Benchmark dose estimation using bootstrap
#' 
#' The function estimates BMD and BMDL using bootstrap based on parametric
#' dose-response models.
#' 
#' BMDL is defined as the 5th percentile in the bootstrap distribution.
#' 
#' Bootstrapping with the argument boot = "nonparametric" is done by sampling
#' with replacement from the original data set. Bootstrapping with the argument
#' boot = "parametric" is done by sampling from norm(mean(Y_i),sd(Y_0)),
#' assuming equal variance between groups, in case of continuous data. For
#' binomial data, each bootstrap data set is sampled from binom(N_i,Y_i/N_i).
#' In case of Y_i = 0 or Y_i = N_i shrinkage is used to avoid that the
#' resampling always produces 0 or 1, respectively. In this case data is
#' sampled from binom(N_i,(Y_i+1/4)/(N_i+1/2)). Bootstrapping with argument
#' boot = "semiparametric" is done by sampling with replacement from the
#' residuals.
#' 
#' All sampling is made within dose groups.
#' 
#' @param object object of class \code{drc}
#' @param bmr numeric value of benchmark response level for which to calculate
#' the benchmark dose
#' @param R number of bootstrap samples. Default is 1000
#' @param bootType character string specifying type of bootstrap used.
#' "nonparametric" (default), "semiparametric" or "parametric". "Semiparametric
#' is only available for continuous data and "nonparametric" is at present the
#' only option for count data. See details below
#' @param bmdType Type of estimate for BMD. Default is "orig" the bmd estimate
#' from the original data set. Other choices are "mean" - the mean of the
#' bootstrap estimates, or "median" - the median of the bootstrap estimates
#' @param backgType character string specifying how the background level is
#' specified. For binomial and count data the options are "modelBased" and
#' "absolute". For continuous data the options are "modelBased", "absolute",
#' "hybridSD" and "hybridPercentile"
#' 
#' "modelBased" - the background level is obtained from the model as the level
#' for dose 0: p0 = f(0)
#' 
#' "absolute" - the background level is specified by the user through the backg
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
#' for "absolute" background level for binomial response (1 for decreasing
#' dose-response models), 2 SD for "hybridSD" background and 0.9 for
#' "hybridpercentile"
#' @param controlSD numeric value specifying the standard deviation of the
#' control group (used in the hybrid approach). If not specified the SD for the
#' control group will be estimated as the square root of the residual variance
#' (assuming variance homogeneity).
#' @param def character string specifying the definition of the benchmark dose
#' to use in the calculations. "excess" , "additional" and "point" are for
#' binomial response whereas "relative", "extra", "added", "hybridExc" (excess
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
#' continuous response
#' 
#' "relative" - BMR is defined as: BMR = (f(BMD) - p0)/p0.  Works for count and
#' continuous response
#' 
#' "extra" - BMR is defined as: BMR = (f(BMD) - p0)/(f(Inf) - p0).  Works for
#' count and continuous response
#' 
#' "added" - BMR is defined as: BMR= f(BMD) + p0.  Works for count and
#' continuous response
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
#' @param bootInterval character string indicating type of bootstrap interval
#' used for estimating BMDL. "Percentile" (default) or "BCa" (Bias corrected
#' and adjusted)
#' @param display logical. If TRUE the results are displayed; otherwise they
#' are not
#' @param level numeric value specifying the levle of the confidence interval
#' underlying BMDL. Default is 0.95
#' @return A list of three elements: Results contain the estimated BMD and
#' BMDL, bootEst is a vector af all of the estimated BMD values from each
#' bootstrap sample, Interval gives BMDL and BMDU, which is identical to the
#' confidence interval for the percentile interval approach.
#' @author Signe M. Jensen
#' @keywords bootstrap
#' @examples
#' 
#' ## Data on root length in ryegrass after exposure to ferulic acid
#' require(drc)
#' require(drcData)
#' data(ryegrass)
#' 
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#' 
#' ## BMD using the hybrid method, background risk is 2 SD, hybrid definition using excess risk
#' bmdBoot(ryegrass.m1, 0.05, backgType = "hybridSD", def = "hybridAdd", R = 50)
#' 
#' ## BMD from the same definitions but using parametric bootstrap
#' bmdBoot(ryegrass.m1, 0.05, backgType = "hybridSD", def = "hybridAdd", bootType="parametric",R = 50)
#' 
#' @export
bmdBoot <- function(object, bmr, R=1000, bootType="nonparametric", bmdType = "orig",
                    backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
                    backg=NA, 
                    controlSD=NA,
                    def = c("excess", "additional", 
                            "relative", "extra", "added", "hybridExc", "hybridAdd", "point"),
                    respTrans = c("none", "log", "sqrt"),
                    bootInterval = c("percentile","BCa"),
                    display=TRUE, level=0.95){
  # Assertions
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
    stop(paste("backgType is missing", sep=""))
  }
  if (!(def %in% c("excess", "additional", "relative", "extra", "added", "hybridExc", "hybridAdd", "point"))) {
    stop(paste("Could not recognize def", sep=""))
  }
  if (!(backgType %in% c("modelBased","absolute","hybridSD","hybridPercentile"))) {
    stop(paste("Could not recognize backgType", sep=""))
  }
  
  if (identical(object$type,"binomial") & bootType=="semiparametric") {
    stop(paste("\"Semiparametric bootstrap does not work for quantal data", sep=""))
  }
  if (object$type %in% c("Poisson","negbin1","negbin2") & bootType!="nonparametric") {
    stop(paste("\"",object$type,"\" only works with nonparametric bootstrap", sep=""))
  }
  
  respTrans <- match.arg(respTrans)
  
  if(inherits(object$fct, "braincousens") & is.null(object$fct$fixed)){
    if(object$fct$name == "BC.4"){
      object$fct$fixed <- c(NA, 0, NA, NA, NA)
    } else if(object$fct$name == "BC.5"){
      object$fct$fixed <- c(NA, NA, NA, NA, NA)
    }
  }
  
  # Set level
  level <- 1-2*(1-level)
  
  # bmd on original data
  bmd0 <- bmd(object = object, bmr = bmr, backgType = backgType,
              backg=backg, controlSD=controlSD, def = def, respTrans = respTrans,
              interval = "delta", display = FALSE)
  
  get.drm.list <- function(tmp.data){
    if(ncol(object$parmMat) == 1){
      drm.list.tmp <- lapply(tmp.data, function(x){
        try(eval(substitute(drm(formula0, data = x, type = object$type, fct = object[["fct"]],
                                control = drmc(noMessage = TRUE)),
                            list(formula0 = object$call$formula)
                            )), TRUE)
      }
      )
    } else if(is.null(object$call$pmodels)){
      drm.list.tmp <- lapply(tmp.data, function(x){
        # if(object$type != "binomial"){
        #   x[[as.character(object$call$curveid)]] <- x[[paste0("orig.", as.character(object$call$curveid))]]
        # }
        try(
          eval(substitute(drm(object$call$formula, weights = weights0, curveid = curveid0,
                              data = x, type = object$type, fct = object$fct, control = drmc(noMessage = TRUE)),
                          list(weights0 = object$call$weights,
                               curveid0 = object$call$curveid)
          )),
          silent = TRUE)
      })
    } else {
      drm.list.tmp <- lapply(tmp.data, function(x){
        # if(object$type != "binomial"){
        #   x[[as.character(object$call$curveid)]] <- x[[paste0("orig.", as.character(object$call$curveid))]]
        # }
        try(
          eval(substitute(drm(formula0, weights = weights0, curveid = curveid0, pmodels = pmodels0,
                              data = x, type = object$type, fct = object$fct, control = drmc(noMessage = TRUE)),
                          list(formula0 = object$call$formula,
                               weights0 = object$call$weights,
                               curveid0 = object$call$curveid,
                               pmodels0 = object$call$pmodels)
          )),
          silent = TRUE)
      })
    }
  }
  
  get.bmd <- function(x){
    if(ncol(object$parmMat) == 1){
      try(bmd(x, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
              display=FALSE, level=level)[["Results"]][1],TRUE)
    } else {
      try(bmd(x, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
              display=FALSE, level=level)[["Results"]][colnames(object$parmMat), 1],TRUE)
    }
  }
 
  if (object$type %in% c("binomial","continuous")) {
    
    tmp.data <- bootDataGen(object,R,bootType,aggregated=FALSE)
    
    drm.list.tmp <- get.drm.list(tmp.data)
    list.condition <- sapply(drm.list.tmp, function(x) class(x)=="drc")
    drm.list  <- drm.list.tmp[list.condition]
    
    bmd.list.try <- lapply(drm.list,get.bmd)
    bmd.list <- suppressWarnings(bmd.list.try[!sapply(bmd.list.try, function(x) any(is.na(as.numeric(x))))])
  }
  
  if (object$type %in% c("Poisson","negbin1","negbin2")) {
    tmp.data <- bootDataGen(object,R,bootType,aggregated=FALSE)
  
    # drm.list.tmp <- lapply(tmp.data, function(x){
    #   try(drm(object$call$formula, data = x, type = object$type, weights=weights, fct = object[["fct"]]),TRUE)}
    # )
    drm.list.tmp <- get.drm.list(tmp.data)
    list.condition <- sapply(drm.list.tmp, function(x) class(x)=="drc")
    drm.list  <- drm.list.tmp[list.condition]
    
    # bmd.list <- lapply(drm.list,function(x){
    #   bmd(x, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
    #       display=FALSE, level=level)[["Results"]][1]}
    # )
    bmd.list.try <- lapply(drm.list,get.bmd)
    bmd.list <- bmd.list.try[!sapply(bmd.list.try, function(x) any(is.na(x)))]
  }

  
  if(identical(object$type, "continuous")){
    if(identical(bootInterval,"BCa")){
      jackData <- list()
      for(i in 1:(dim(object$data)[1])){
        jackData[[i]] <- object$origData[-i,]
      }
      # bootJack.drm.tmp <- lapply(jackData, function(x){
      #   try(drm(object$call$formula, data = x, fct = object[["fct"]]),TRUE)
      # })
      bootJack.drm.tmp <- get.drm.list(jackData)
      list.condition <- sapply(bootJack.drm.tmp, function(x) class(x)=="drc")
      bootJack.drm<- bootJack.drm.tmp[list.condition]
      
      # bootJack <- sapply(bootJack.drm, function(x){
      #   bmd(x, bmr, backgType = backgType, backg = backg, def = def, controlSD=controlSD, respTrans = respTrans,
      #       interval = "delta", display=FALSE, level=level)$Results[1]
      # }
      # )
      
      bootJack.list.try <- lapply(bootJack.drm, get.bmd)
      bootJack.list <- bootJack.list.try[!sapply(bootJack.list.try, function(x) any(is.na(x)))]
      
      # use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
      #                display=FALSE, level=level)[["Results"]][1]
      use.bmd <- get.bmd(object)
      if(ncol(object$parmMat)==1){
        BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = object$data, unlist(bmd.list), unlist(bootJack.list), level=level)[1])
      } else {
        BCaBMDL <- sapply(1:ncol(object$parmMat), function(i) as.numeric(BCa(obs = use.bmd[i], data = object$data, 
                                                                             sapply(bmd.list, function(x) x[i]), 
                                                                             sapply(bootJack.list, function(x) x[i]), level = level)[1]))
      }
      
    }
  }
  if(identical(object$type, "binomial")){
    if(identical(bootInterval, "BCa")){
      # data.str <- object$data
      # data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
      # data.e<-expandBinomial(data.str, 
      #                        number = "number",
      #                        total = "weights",
      #                        dose = as.character(object$call$formula[[3]]))
      
      if(ncol(object$parmMat) == 1){
        data.str <- object$data
        data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
        data.e<-expandBinomial(data.str, 
                               number = "number",
                               total = "weights",
                               dose = as.character(object$call$formula[[3]]))
        df <- data.frame(data.e[,as.character(object$call$formula[[3]])],
                         data.e[,"number"],
                         data.e[,"weights"])
        colnames(df) <- c(as.character(object$call$formula[[3]]),
                          as.character(object$call$formula[[2]])[[2]],
                          as.character(object$call$formula[[2]])[[3]])
        df$rowNum <-1:nrow(df)
      } else {
        data.str <- object$data
        data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
        data.e<-expandBinomial(data.str, 
                               number = "number",
                               total = "weights",
                               dose = as.character(object$call$formula[[3]]),
                               curveid = as.character(object$call$curveid))
        df <- data.frame(data.e[,as.character(object$call$formula[[3]])],
                         data.e[,"number"],
                         data.e[,"weights"],
                         data.e[,as.character(object$call$curveid)])
        colnames(df) <- c(as.character(object$call$formula[[3]]),
                          as.character(object$call$formula[[2]])[[2]],
                          as.character(object$call$formula[[2]])[[3]],
                          as.character(object$call$curveid))
      }
      
      jackData <- list()
      for(i in 1:(dim(data.e)[1])){
        # jackData[[i]] <- data.e[-i,]
        jackData[[i]] <- df[-i,]
      }
      
      # bootJack.drm.tmp <- lapply(jackData, function(x){
      #   try(drm(as.formula(paste0("number~", as.character(object$call$formula[[3]]))), data = x, type = "binomial", fct = object[["fct"]]),TRUE) # number~dose
      # })
      bootJack.drm.tmp <- get.drm.list(jackData)
      list.condition <- sapply(bootJack.drm.tmp, function(x) class(x)=="drc")
      bootJack.drm<- bootJack.drm.tmp[list.condition]
      
      # bootJack <- sapply(bootJack.drm, function(x){
      #   bmd(x, bmr, backgType = backgType, backg = backg, def = def, controlSD=controlSD, respTrans = respTrans,
      #       interval = "delta", display=FALSE, level=level)$Results[1]
      # }
      # )
      
      bootJack.list.try <- lapply(bootJack.drm, get.bmd)
      bootJack.list <- bootJack.list.try[!sapply(bootJack.list.try, function(x) any(is.na(x)))]
      
      # use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
      #                display=FALSE, level=level)[["Results"]][1]
      # BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = data.e, unlist(bmd.list), bootJack, level=level)[1])
      
      use.bmd <- get.bmd(object)
      if(ncol(object$parmMat)==1){
        BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = df, unlist(bmd.list), unlist(bootJack.list), level=level)[1])
      } else {
        BCaBMDL <- sapply(1:ncol(object$parmMat), function(i) as.numeric(BCa(obs = use.bmd[i], data = df, 
                                                                             sapply(bmd.list, function(x) x[i]), 
                                                                             sapply(bootJack.list, function(x) x[i]), level = level)[1]))
      }
      
    }
  }
  if(object$type %in% c("Poisson","negbin1","negbin2")){
    if(identical(bootInterval,"BCa")){
      jackData <- list()
      for(i in 1:(dim(object$data)[1])){
        jackData[[i]] <- object$data[-i,]
      }
      # bootJack.drm.tmp <- lapply(jackData, function(x){
      #   try(drm(object$call$formula, data = x, type = object$type, weights=weights, fct = object[["fct"]]),TRUE)
      # })
      bootJack.drm.tmp <- get.drm.list(jackData)
      list.condition <- sapply(bootJack.drm.tmp, function(x) class(x)=="drc")
      bootJack.drm<- bootJack.drm.tmp[list.condition]
      
      # bootJack <- sapply(bootJack.drm, function(x){
      #   bmd(x, bmr, backgType = backgType, backg = backg, def = def, controlSD=controlSD, interval = "delta", respTrans = respTrans,
      #       display=FALSE, level=level)$Results[1]
      # }
      # )
      
      bootJack.list.try <- lapply(bootJack.drm, get.bmd)
      bootJack.list <- bootJack.list.try[!sapply(bootJack.list.try, function(x) any(is.na(x)))]
      
      # use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
      #                display=FALSE, level=level)[["Results"]][1]
      # BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = object$data, bootSample=unlist(bmd.list), 
      #                           bootjack=bootJack, level=level)[1])
      use.bmd <- get.bmd(object)
      if(ncol(object$parmMat)==1){
        BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = object$data, unlist(bmd.list), unlist(bootJack.list), level=level)[1])
      } else {
        BCaBMDL <- sapply(1:ncol(object$parmMat), function(i) as.numeric(BCa(obs = use.bmd[i], data = object$data, 
                                                                             sapply(bmd.list, function(x) x[i]), 
                                                                             sapply(bootJack.list, function(x) x[i]), level = level)[1]))
      }
    }
  }
  if(bmdType == "orig"){
    if(ncol(object$parmMat) == 1){
      use.bmd <- bmd0[["Results"]][1]
    } else {
      use.bmd <- bmd0[["Results"]][colnames(object$parmMat), 1]
    }
  } else if(bmdType == "mean"){
    if(ncol(object$parmMat) == 1){
      use.bmd <- mean(unlist(bmd.list))
    } else {
      use.bmd <- apply(matrix(unlist(bmd.list), ncol = ncol(object$parmMat), byrow = TRUE), 2, mean, na.rm = TRUE)
    }
  } else if(bmdType == "median"){
    if(ncol(object$parmMat) == 1){
      use.bmd <- quantile(unlist(bmd.list),0.5)  
    } else {
      use.bmd <- apply(matrix(unlist(bmd.list), ncol = ncol(object$parmMat), byrow = TRUE), 2, quantile, p= 0.5, na.rm = TRUE)
    }
  } 
  
  if(ncol(object$parmMat) == 1){
    BMDL <- ifelse(identical(bootInterval,"BCa"), BCaBMDL, quantile(unlist(bmd.list),c(1-level),na.rm = TRUE))
    BMDU <- ifelse(identical(bootInterval,"BCa"), "Not available for BCa bootstrap", quantile(unlist(bmd.list),c(level),na.rm = TRUE))
  } else {
    if(identical(bootInterval, "BCa")){
      BMDL <- BCaBMDL
      BMDU <- rep("Not available for BCa bootstrap", ncol(object$parmMat))
    } else {
      BMDL <- apply(matrix(unlist(bmd.list), ncol = ncol(object$parmMat), byrow = TRUE), 
                    2, quantile, p= c(1-level), na.rm = TRUE)
      BMDU <- apply(matrix(unlist(bmd.list), ncol = ncol(object$parmMat), byrow = TRUE), 2, quantile, p= c(level), na.rm = TRUE)
    }
  }
  
  resMat <- matrix(NA,ncol(object$parmMat),2)
  resMat[,1] <- use.bmd
  resMat[,2] <- BMDL
  colnames(resMat) <- c("BMD", "BMDL")
  
  intMat <- matrix(NA,ncol(object$parmMat),2)
  intMat[,1] <- BMDL
  intMat[,2] <- BMDU
  colnames(intMat) <- c("Lower", "Upper")
  
  if(ncol(object$parmMat) == 1){
    rownames(resMat) <- c("")
    rownames(intMat) <- c("")
    bootEst = unlist(bmd.list)
    used.Boot <- sum(!is.na(bootEst))
  } else {
    rownames(resMat) <- colnames(object$parmMat)
    rownames(intMat) <- colnames(object$parmMat)
    bootEst <- matrix(unlist(bmd.list), byrow = TRUE, ncol = ncol(object$parmMat), dimnames = list(NULL, colnames(object$parmMat)))
    used.Boot <- sum(sapply(bmd.list, function(x) all(!is.na(x))))
  }
  
  if(display){
    print(resMat)
  }
  
  used.Boot <- sum(sapply(bmd.list, function(x) all(!is.na(x))))
  
  resBMD<-list(Results = resMat,
               Boot.samples.used = used.Boot,
               bootEst = bootEst,
               interval = intMat)
  class(resBMD) <- "bmd"
  invisible(resBMD)
}
