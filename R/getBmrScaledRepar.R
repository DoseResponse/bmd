getBmrScaledRepar <- function(cdVec, slope, bmr, backgType, 
                             backg=NA, #controlSD=NA,
                             def,
                             respTrans, respType){
  
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
    cat('respTrans "', respTrans, '" not defined. Options are: "none", "log", and "sqrt".')
    break
  }
  
  # useSD
  # if(def %in% c("hybridAdd","hybridExc")){
  #   useSD <- ifelse(!is.na(controlSD),controlSD,sqrt(summary(object)$resVar))
  # }
  
  # INCREASING CURVE
  if(slope == "increasing"){
    # BACKGROUND
    if (identical(backgType,"modelBased")) {
      background <- cdVec[1] # c parameter
    } else if (identical(backgType,"absolute") & !(def %in% c("relative","extra"))) {
      background <- backg
    } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
      background <- ifelse(is.na(backg),0,backg)
    } else {
      background <- ifelse(is.na(backg),1-0.9,1-backg)
    }
    # if (identical(backgType,"hybridSD")) {
    #   background <- ifelse(is.na(backg), 1-pnorm(2), 1-pnorm(backg))
    # } 
    # if (identical(backgType,"absolute") & 
    #     (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
    #   if(!identical(respTrans, "none") & !is.na(backg)){
    #     backg <- h(backg)
    #   }
    #   background <- ifelse(is.na(backg), 
    #                        1-pnorm(2),
    #                        1-pnorm((backg-cdVec[1])/useSD))
    # } 
    
    # Apply inverse transformation to calculate bmrScaled on original scale
    if(!identical(respTrans, "none") & !(def %in% c("hybridExc","hybridAdd"))){
      background <- hInv(background)
    }
    
    # BMRSCALED
    if (identical(respType, "binomial")) {
      bmrScaled <- switch(def, 
                          excess = bmr * (1 - background) + background, 
                          additional = bmr + background, 
                          point = bmr)
    }
    if (identical(respType, "binomial") & (def %in% c("relative","added","extra", "hybridExc","hybridAdd"))) {
      stop(paste("\"",def, "\" is not available for quantal data", sep=""))
    }
    if (respType %in% c("Poisson","negbin1","negbin2")) {
      bmrScaled <- switch(def, 
                          relative = bmr * background + background,
                          added = bmr + background,
                          extra = bmr*diff(cdVec) + background,
                          point = bmr)
    }
    if (respType %in% c("Poisson","negbin1","negbin2") & (def %in% c("excess","additional","hybridExc","hybridAdd"))) {
      stop(paste("\"",def, "\" is not available for count data", sep=""))
    }
    # if (def %in% c("excess","additional") & (backgType %in% c("hybridSD", "hybridPercentile"))) {
    #   stop(paste("\"",def, "\" does not work with backgType = \"", backgType,"\"", sep=""))
    # }
    if (identical(respType, "continuous")) {
      bmrScaled <- switch(def, 
                          relative = bmr * background + background, 
                          added = bmr + background,
                          point = bmr,
                          extra = bmr*diff(cdVec) + background#,
                          # hybridAdd = useSD * 
                          #   (qnorm(1 - background) - qnorm(1 - (background + bmr))) + 
                          #   cdVec[1],
                          # hybridExc = useSD * 
                          #   (qnorm(1 - background) - qnorm(1 - background - (1 - background)*bmr)) + 
                          #   cdVec[1]
                          )
    }
    
    if (identical(respType, "continuous") & (def %in% c("excess", "additional"))) {
      stop(paste("\"",def, "\" is not available for continuous data", sep=""))
    }
    
    # Transform back to transformed scale, so bmrScaled is same scale as model
    if(!identical(respTrans, "none") & !(def %in% c("hybridExc","hybridAdd"))){
      bmrScaled <- h(bmrScaled)
    }
    
  } 
  # DECREASING CURVE
  else if(slope == "decreasing"){
    # BACKGROUND
    if (identical(backgType,"modelBased")) {
      background <- cdVec[2] # d parameter
    } else if (identical(backgType,"absolute") & !(def %in% c("relative","extra"))) {
      background <- backg
    } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
      background <- ifelse(is.na(backg),1,backg)
    } else {
      background <- ifelse(is.na(backg),0.1,backg)
    }
    # if (identical(backgType,"hybridSD")) {
    #   background <- ifelse(is.na(backg), pnorm(-2), pnorm(-backg))
    # } 
    # if (identical(backgType,"absolute") & 
    #     (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
    #   if(!identical(respTrans, "none") & !is.na(backg)){
    #     backg <- h(backg)
    #   }
    #   background <- ifelse(is.na(backg), 
    #                        pnorm(-2),
    #                        pnorm((backg-cdVec[2])/useSD))
    # } 
    
    # Apply inverse transformation to calculate bmrScaled on original scale
    if(!identical(respTrans, "none") & !(def %in% c("hybridExc","hybridAdd"))){
      background <- hInv(background)
    }
    
    # BMRSCALED
    if (identical(respType, "binomial")) {
      bmrScaled <- switch(def, 
                          excess = background - bmr * background, 
                          additional = background - bmr, 
                          point = bmr)
    }
    if (identical(respType, "binomial") & (def %in% c("relative","added","extra", "hybridExc","hybridAdd"))) {
      stop(paste("\"",def, "\" is not available for quantal data", sep=""))
    }
    if (respType %in% c("Poisson","negbin1","negbin2")) {
      bmrScaled <- switch(def, 
                          relative = background - bmr * background,
                          added = background - bmr,
                          extra = background - bmr*diff(cdVec),
                          point = bmr)
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
                          extra = background - bmr*diff(cdVec)#,
                          # hybridAdd = useSD * 
                          #   (qnorm(background) - qnorm(background + bmr)) + 
                          #   cdVec[2],
                          # hybridExc = useSD * 
                          #   (qnorm(background) - qnorm(bmr + (1-bmr) * background)) + 
                          #   cdVec[2]
                          )
    }
    
    if (identical(respType, "continuous") & (def %in% c("excess", "additional"))) {
      stop(paste("\"",def, "\" is not available for continuous data", sep=""))
    }
    
    # Transform back to transformed scale, so bmrScaled is same scale as model
    if(!identical(respTrans, "none") & !(def %in% c("hybridExc","hybridAdd"))){
      bmrScaled <- h(bmrScaled)
    }
  }
  
  bmrScaled
}