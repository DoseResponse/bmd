bmd<-function (object, bmr, backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
               backg=NA, controlSD=NA,
               def = c("excess", "additional", 
                       "relative", "extra", "added", "hybridExc", "hybridAdd", "point"), 
              interval = "delta", sandwich.vcov=FALSE, display = TRUE, level=0.95) 
{
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
  slope <- ifelse(as.numeric(predict(object,data.frame(0))-predict(object,data.frame(Inf)))>0,"decreasing","increasing")
    if( identical(slope,"increasing" )) {
    f0 <- ifelse(!is.na(coef(object)["c:(Intercept)"]), coef(object)["c:(Intercept)"], predict(object,data.frame(0)))
    } else {
    f0 <- ifelse(!is.na(coef(object)["d:(Intercept)"]), coef(object)["d:(Intercept)"], predict(object,data.frame(0)))
    }
  
  if(def %in% c("hybridAdd","hybridExc")){
  useSD <- ifelse(!is.na(controlSD),controlSD,sqrt(summary(object)$resVar))
  }
  
  if(identical(slope,"increasing" )) {
      if (identical(backgType,"modelBased")) {
        background <- f0
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
      background <- ifelse(is.na(backg), 
                           1-pnorm(2),
                           1-pnorm((backg-f0)/useSD))
    } 
    } 
  if(identical(slope,"decreasing" )) {
      if (identical(backgType,"modelBased")) {
        background <- f0
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
      background <- ifelse(is.na(backg), 
                           pnorm(-2),
                           pnorm((backg-f0)/useSD))
    }
    }
    
    def <- match.arg(def)
    respType <- object$type
    if (identical(respType, "binomial")) {
      if(identical(slope,"increasing" )) {
        bmrScaled <- switch(def, 
                            excess = bmr * (1 - background) + background, 
                            additional = bmr + background, 
                            point = bmr)
      } else {
        bmrScaled <- switch(def, 
                            excess = background - bmr * background, 
                            additional = background - bmr, 
                            point = bmr)
      } 
      typeVal <- "absolute"
    }
    if (identical(respType, "binomial") & (def %in% c("relative","added","extra", "hybridExc","hybridAdd"))) {
      stop(paste("\"",def, "\" is not available for quantal data", sep=""))
    }
    if (respType %in% c("Poisson","negbin1","negbin2")) {
      if(identical(slope,"increasing" )) {
        bmrScaled <- switch(def, 
                            relative = bmr * background + background,
                            added = bmr + background,
                            extra = bmr*abs(diff(predict(object, data.frame(c(0, Inf))))) + background,
                            point = bmr)
      } else {
        bmrScaled <- switch(def, 
                            relative = background - bmr * background,
                            added = background - bmr,
                            extra = background - bmr*abs(diff(predict(object, data.frame(c(0, Inf))))),
                            point = bmr)
      } 
      typeVal <- "absolute"
    }
    if (respType %in% c("Poisson","negbin1","negbin2") & (def %in% c("excess","additional","hybridExc","hybridAdd"))) {
      stop(paste("\"",def, "\" is not available for count data", sep=""))
    }
    if (def %in% c("excess","additional") & (backgType %in% c("hybridSD", "hybridPercentile"))) {
      stop(paste("\"",def, "\" does not work with backgType = \"", backgType,"\"", sep=""))
    }
    if (identical(respType, "continuous")) {
      if(identical(slope,"increasing" )) {
        bmrScaled0 <- switch(def, 
                             relative = bmr * background + background, 
                             added = bmr + background,
                             point = bmr,
                             extra = bmr*abs(diff(predict(object, data.frame(c(0, Inf)))))
                                     + background,
                             hybridAdd = useSD * 
                                   (qnorm(1 - background) - qnorm(1 - (background + bmr))) + 
                                   f0,
                             hybridExc = useSD * 
                                  (qnorm(1 - background) - qnorm(1 - background - (1 - background)*bmr)) + 
                                   f0)
      } else {
        bmrScaled0 <- switch(def, 
                             relative = background - bmr * background, 
                             added = background - bmr,
                             point = bmr,
                             extra = background - bmr*abs(diff(predict(object, data.frame(c(0, Inf))))),
                             hybridAdd = useSD * 
                               (qnorm(background) - qnorm(background + bmr)) + 
                               f0,
                             hybridExc = useSD * 
                               (qnorm(background) - qnorm(bmr - background*bmr + background)) + 
                               f0)
      }
        bmrScaled <- as.numeric(format(bmrScaled0, digits = 5))
        typeVal <- "absolute"
    }
    
    if (identical(respType, "continuous") & (def %in% c("excess", "additional"))) {
      stop(paste("\"",def, "\" is not available for continuous data", sep=""))
    }
    
    if(interval=="delta"){
    if(sandwich.vcov==TRUE){
    resMat <- ED(object, bmrScaled, interval = "delta", 
                level = 1-2*(1-level), type = typeVal, vcov. = sandwich, display = FALSE)[, 
                  c("Estimate", "Lower"), drop = FALSE]
    }
    if(sandwich.vcov==FALSE){
    resMat <- ED(object, bmrScaled, interval = interval, 
                level = 1-2*(1-level), type = typeVal, vcov. = vcov, display = FALSE)[, 
                  c("Estimate", "Lower"), drop = FALSE]
    }
    colnames(resMat) <- c("BMD", "BMDL")
    rownames(resMat) <- c("")
    if(sandwich.vcov==TRUE){
      bmdInterval <- ED(object, bmrScaled, interval = interval, 
                        level = 1-2*(1-level), type = typeVal, vcov. = sandwich, display = FALSE)[, 
                                                                                     c("Lower", "Upper"), drop = FALSE]
    }
    if(sandwich.vcov==FALSE){
      bmdInterval <- ED(object, bmrScaled, interval = interval, 
                        level = 1-2*(1-level), type = typeVal, vcov. = vcov, display = FALSE)[, 
                                                                                        c("Lower", "Upper"), drop = FALSE]
    }
    colnames(bmdInterval) <- c("Lower CI", "Upper CI")
    rownames(bmdInterval) <- c("")
    }
    
    if(interval=="inv"){
      resMat <- inv.bmd(object, bmr, level = 1-2*(1-level), slope=slope, backgType=backgType,
                             backg=backg, catLev=NA, extFactor=10, def=def, useSD=useSD, sandwich.vcov=sandwich.vcov)[, 
                                                                                               c("BMD","BMDL"), drop = FALSE]
      colnames(resMat) <- c("BMD", "BMDL")
      rownames(resMat) <- c("")
      bmdInterval <- inv.bmd(object, bmr, level = 1-2*(1-level), slope=slope, backgType=backgType,
                                  backg=backg, catLev=NA, extFactor=10, def=def, useSD=useSD, sandwich.vcov=sandwich.vcov)[, 
                                                                                                    c("BMDL", "BMDU"), drop = FALSE]
      colnames(bmdInterval) <- c("Lower CI", "Upper CI")
      rownames(bmdInterval) <- c("")
    }
    bmdSE <- matrix(NA,1,1)
    if(sandwich.vcov==TRUE){
      bmdSE[1,1] <- ifelse(identical(interval,"delta"),
                           ED(object, bmrScaled, interval = "delta", 
                              level = 1-2*(1-level), type = typeVal, vcov. = sandwich, 
                              display = FALSE)[, c("Std. Error"), drop = FALSE],
                           NA)
    }
    if(sandwich.vcov==FALSE){
    bmdSE[1,1] <- ifelse(identical(interval,"delta"),
                         ED(object, bmrScaled, interval = "delta", 
                            level = 1-2*(1-level), type = typeVal, vcov. = vcov, 
                            display = FALSE)[, c("Std. Error"), drop = FALSE],
                         NA)
    }
    colnames(bmdSE) <- c("Std. Error")
    rownames(bmdSE) <- c("")
    if (display) {
      print(resMat)
    }
    
    resBMD<-list(Results = resMat,
                 bmrScaled = bmrScaled,
                 interval = bmdInterval,
                 SE = bmdSE)
    class(resBMD) <- "bmd"
    invisible(resBMD) 
    }



