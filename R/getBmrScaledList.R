getBmrScaledList <- function(object, bmr, backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"), 
                             backg=NA, controlSD=NA,
                             def = c("excess", "additional", 
                                     "relative", "extra", "added", "hybridExc", "hybridAdd", "point"),
                             respTrans = c("none", "log", "sqrt")){
  curveLevels <- colnames(object$parmMat)
  nCurves <- ncol(object$parmMat)
  
  bmrScaledMat <- matrix(0, nrow = nCurves, ncol = 1, dimnames = list(curveLevels, NULL))
  dBmrScaledMat <- matrix(0, nrow = nrow(object$parmMat), ncol = ncol(object$parmMat), dimnames = dimnames(object$parmMat))
  
  backgType <- match.arg(backgType)
  def <- match.arg(def)
  respType <- object$type
  
  if(identical(respTrans, "log")){
    h <- log
    dh <- function(x) 1/x
    hInv <- exp
    dhInv <- exp
  } else if(identical(respTrans, "sqrt")){
    h <- sqrt
    dh <- function(x) 1/(2*sqrt(x))
    hInv <- function(x) x^2
    dhInv <- function(x) 2*x
  } else if(!identical(respTrans, "none")){
    stop(paste0('respTrans "', respTrans, '" not defined. Options are: "none", "log", and "sqrt".'))
  }
  
  zero_one_vector <- function(onePos, vecLength){
    vec <- numeric(vecLength)
    vec[onePos] <- 1
    vec
  }
  
  for(iCurve in 1:nCurves){
    # Vectors of parameters
    fullParmVec <- object$fct$fixed
    notFixed <- is.na(object$fct$fixed)
    fullParmVec[notFixed] <- object$parmMat[,iCurve]
    dFullParmVec <- numeric(length(fullParmVec))
    
    # SLOPE
    curveMatColIndex <- unique(object$dataList$curveid) == colnames(object$parmMat)[iCurve]
    slope <- drop(ifelse(object$curve[[1]](0)[,curveMatColIndex]-object$curve[[1]](Inf)[,curveMatColIndex]>0,"decreasing","increasing"))
    if(is.na(object$curve[[1]](0)[,curveMatColIndex]-object$curve[[1]](Inf)[,curveMatColIndex])){
      slope <- drop(ifelse(object$curve[[1]](0.00000001)[,curveMatColIndex]-object$curve[[1]](100000000)[,curveMatColIndex]>0,"decreasing","increasing"))
    }
    
    # useSD
    if(def %in% c("hybridAdd","hybridExc")){
      useSD <- ifelse(!is.na(controlSD),controlSD,sqrt(summary(object)$resVar))
    }
    
    # INCREASING CURVE
    if(slope == "increasing"){
      # BACKGROUND
      dBackground <- numeric(length(fullParmVec))
      f0_par_ind <- which.min(fullParmVec[2:3]) + 1
      if (identical(backgType,"modelBased")) {
        background <- fullParmVec[f0_par_ind]
        dBackground[f0_par_ind] <- 1
      } else if (identical(backgType,"absolute") & !(def %in% c("relative","extra"))) {
        background <- backg
      } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
        background <- ifelse(is.na(backg),0,backg)
      } else {
        background <- ifelse(is.na(backg),1-0.9,1-backg)
      }
      if (identical(backgType,"hybridSD")) {
        background <- ifelse(is.na(backg), 1-pnorm(2), 1-pnorm(backg))
      } 
      if (identical(backgType,"absolute") & 
          (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
        if(!identical(respTrans, "none") & !is.na(backg)){
          backg <- h(backg)
        }
        background <- ifelse(is.na(backg), 
                             1-pnorm(2),
                             1-pnorm((backg-fullParmVec[f0_par_ind])/useSD))
        dBackground[f0_par_ind] <- ifelse(is.na(backg), 
                             0,
                             dnorm((backg-fullParmVec[f0_par_ind])/useSD)/useSD)
      } 

      # Apply inverse transformation to calculate bmrScaled on original scale
      if(!identical(respTrans, "none") & !(def %in% c("hybridExc","hybridAdd"))){
        dBackground <- dhInv(background) * dBackground
        background <- hInv(background)
      }
      
      # BMRSCALED
      if (identical(respType, "binomial")) {
        bmrScaled <- switch(def, 
                            excess = bmr * (1 - background) + background, 
                            additional = bmr + background, 
                            point = bmr)
        dFullParmVec <- switch(def, 
                            excess = bmr * ( - dBackground) + dBackground, 
                            additional = dBackground, 
                            point = rep(0, length(fullParmVec)))
      }
      if (identical(respType, "binomial") & (def %in% c("relative","added","extra", "hybridExc","hybridAdd"))) {
        stop(paste("\"",def, "\" is not available for quantal data", sep=""))
      }
      if (respType %in% c("Poisson","negbin1","negbin2")) {
        bmrScaled <- switch(def, 
                            relative = bmr * background + background,
                            added = bmr + background,
                            extra = bmr*abs(diff(fullParmVec[2:3])) + background,
                            point = bmr)
        dFullParmVec <- switch(def, 
                            relative = bmr * dBackground + dBackground,
                            added = dBackground,
                            extra = bmr*ifelse(rep(f0_par_ind==2, length(fullParmVec)), c(0,-1,1, rep(0,length(fullParmVec)-3)), c(0,1,-1, rep(0,length(fullParmVec)-3))) + dBackground,
                            point = rep(0, length(fullParmVec)))
      }
      if (respType %in% c("Poisson","negbin1","negbin2") & (def %in% c("excess","additional","hybridExc","hybridAdd"))) {
        stop(paste("\"",def, "\" is not available for count data", sep=""))
      }
      if (def %in% c("excess","additional") & (backgType %in% c("hybridSD", "hybridPercentile"))) {
        stop(paste("\"",def, "\" does not work with backgType = \"", backgType,"\"", sep=""))
      }
      if (identical(respType, "continuous")) {
        bmrScaled <- switch(def, 
                             relative = bmr * background + background, 
                             added = bmr + background,
                             point = bmr,
                             extra = bmr*abs(diff(fullParmVec[2:3])) + background,
                             hybridAdd = useSD * 
                               (qnorm(1 - background) - qnorm(1 - (background + bmr))) + 
                              fullParmVec[f0_par_ind],
                             hybridExc = useSD * 
                               (qnorm(1 - background) - qnorm(1 - background - (1 - background)*bmr)) + 
                              fullParmVec[f0_par_ind])
        dFullParmVec <- switch(def, 
                             relative = bmr * dBackground + dBackground, 
                             added = dBackground,
                             point = rep(0, length(fullParmVec)),
                             extra = bmr*ifelse(rep(f0_par_ind==2, length(fullParmVec)), c(0,-1,1, rep(0,length(fullParmVec)-3)), c(0,1,-1, rep(0,length(fullParmVec)-3))) + dBackground,
                             hybridAdd = useSD * 
                               (-dBackground/dnorm(qnorm(1 - background)) + dBackground/dnorm(qnorm(1 - (background + bmr)))) + 
                               zero_one_vector(f0_par_ind, length(fullParmVec)),
                             hybridExc = useSD * 
                               (-dBackground/dnorm(qnorm(1 - background)) + (1-bmr)*dBackground/dnorm(qnorm(1 - background - (1 - background)*bmr))) + 
                               zero_one_vector(f0_par_ind, length(fullParmVec)))
      }
      
      if (identical(respType, "continuous") & (def %in% c("excess", "additional"))) {
        stop(paste("\"",def, "\" is not available for continuous data", sep=""))
      }
      
      # Transform back to transformed scale, so bmrScaled is same scale as model
      if(!identical(respTrans, "none") & !(def %in% c("hybridExc","hybridAdd"))){
        dFullParmVec <- dh(bmrScaled) * dFullParmVec
        bmrScaled <- h(bmrScaled)
      }
      
    } 
    # DECREASING CURVE
    else if(slope == "decreasing"){
      # BACKGROUND
      dBackground <- numeric(length(fullParmVec))
      f0_par_ind <- which.max(fullParmVec[2:3]) + 1
      
      if (identical(backgType,"modelBased")) {
        background <- fullParmVec[f0_par_ind] # d parameter
        dBackground[f0_par_ind] <- 1
      } else if (identical(backgType,"absolute") & !(def %in% c("relative","extra"))) {
        background <- backg
      } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
        background <- ifelse(is.na(backg),1,backg)
      } else {
        background <- ifelse(is.na(backg),0.1,backg)
      }
      if (identical(backgType,"hybridSD")) {
        background <- ifelse(is.na(backg), pnorm(-2), pnorm(-backg))
      } 
      if (identical(backgType,"absolute") & 
          (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
        if(!identical(respTrans, "none") & !is.na(backg)){
          backg <- h(backg)
        }
        background <- ifelse(is.na(backg), 
                             pnorm(-2),
                             pnorm((backg-fullParmVec[f0_par_ind])/useSD))
        dBackground[f0_par_ind] <- ifelse(is.na(backg), 
                                 0,
                                 -dnorm((backg-fullParmVec[f0_par_ind])/useSD)/useSD)
      }
      
      # Apply inverse transformation to calculate bmrScaled on original scale
      if(!identical(respTrans, "none") & !(def %in% c("hybridExc","hybridAdd"))){
        dBackground <- dhInv(background) * dBackground
        background <- hInv(background)
      }
      
      # BMRSCALED
      if (identical(respType, "binomial")) {
        bmrScaled <- switch(def, 
                            excess = background - bmr * background, 
                            additional = background - bmr, 
                            point = bmr)
        dFullParmVec <- switch(def, 
                               excess = dBackground - bmr * dBackground, 
                               additional = dBackground, 
                               point = rep(0, length(fullParmVec)))
      }
      if (identical(respType, "binomial") & (def %in% c("relative","added","extra", "hybridExc","hybridAdd"))) {
        stop(paste("\"",def, "\" is not available for quantal data", sep=""))
      }
      if (respType %in% c("Poisson","negbin1","negbin2")) {
        bmrScaled <- switch(def, 
                            relative = background - bmr * background,
                            added = background - bmr,
                            extra = background - bmr*abs(diff(fullParmVec[2:3])),
                            point = bmr)
        dFullParmVec <- switch(def, 
                               relative = dBackground - bmr * dBackground,
                               added = dBackground,
                               extra = dBackground - bmr* ifelse(rep(f0_par_ind==2, length(fullParmVec)), c(0,-1,1, rep(0,length(fullParmVec)-3)), c(0,1,-1, rep(0,length(fullParmVec)-3))),
                               point = rep(0, length(fullParmVec)))
      }
      if (respType %in% c("Poisson","negbin1","negbin2") & (def %in% c("excess","additional","hybridExc","hybridAdd"))) {
        stop(paste("\"",def, "\" is not available for count data", sep=""))
      }
      if (def %in% c("excess","additional") & (backgType %in% c("hybridSD", "hybridPercentile"))) {
        stop(paste("\"",def, "\" does not work with backgType = \"", backgType,"\"", sep=""))
      }
      if (identical(respType, "continuous")) {
        bmrScaled <- switch(def, 
                            relative = background - bmr * background, 
                            added = background - bmr,
                            point = bmr,
                            extra = background - bmr*abs(diff(fullParmVec[2:3])),
                            hybridAdd = useSD * 
                              (qnorm(background) - qnorm(background + bmr)) + 
                              fullParmVec[f0_par_ind],
                            hybridExc = useSD * 
                              (qnorm(background) - qnorm(bmr + (1-bmr) * background)) + 
                              fullParmVec[f0_par_ind])
        dFullParmVec <- switch(def, 
                               relative =  dBackground - bmr * dBackground, 
                               added = dBackground,
                               point = rep(0, length(fullParmVec)),
                               extra = dBackground - bmr*ifelse(rep(f0_par_ind==2, length(fullParmVec)), c(0,-1,1, rep(0,length(fullParmVec)-3)), c(0,1,-1, rep(0,length(fullParmVec)-3))),
                               hybridAdd = useSD * 
                                 (dBackground/dnorm(qnorm(background)) - dBackground/dnorm(qnorm(background + bmr))) + 
                                 zero_one_vector(f0_par_ind, length(fullParmVec)),
                               hybridExc = useSD * 
                                 (dBackground/dnorm(qnorm(background)) - (1-bmr)*dBackground/dnorm(qnorm(bmr + (1-bmr) * background))) + 
                                 zero_one_vector(f0_par_ind, length(fullParmVec)))
      }
      
      if (identical(respType, "continuous") & (def %in% c("excess", "additional"))) {
        stop(paste("\"",def, "\" is not available for continuous data", sep=""))
      }
      
      # Transform back to transformed scale, so bmrScaled is same scale as model
      if(!identical(respTrans, "none") & !(def %in% c("hybridExc","hybridAdd"))){
        dFullParmVec <- dh(bmrScaled) * dFullParmVec
        bmrScaled <- h(bmrScaled)
      }
    }
  
    
    # Set values
    bmrScaledMat[iCurve,] <- bmrScaled
    dBmrScaledMat[,iCurve] <- dFullParmVec[notFixed]
  }
  
  list(bmrScaledMat = bmrScaledMat, dBmrScaledMat = dBmrScaledMat)
}
