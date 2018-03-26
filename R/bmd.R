bmd<-function (object, bmr, backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
               backg=NA, 
               def = c("excess", "additional", 
                       "relative", "extra", "added", "hybridExc", "hybridAdd", "point"), 
              interval = "delta", display = FALSE) 
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
    if(identical(slope,"increasing" )) {
      if (identical(backgType,"modelBased")) {
        background <- f0
        } else if (identical(backgType,"absolute") & !(def %in% c("relative","extra"))) {
        background <- backg
        } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
          background <- ifelse(is.na(backg),0,backg)
        } else if (identical(backgType,"hybridSD")) {
        background <- ifelse(is.na(backg), 1-pnorm(2), 1-pnorm(backg))
        } else if (identical(backgType,"absolute") & 
               (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
        background <- ifelse(is.na(backg), 
                             1-pnorm(2),
                             1-pnorm((backg-f0)/sqrt(summary(object)$resVar)))
        } else {
        background <- ifelse(is.na(backg),1-0.9,1-backg)}
    } else {
      if (identical(backgType,"modelBased")) {
        background <- f0
        } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
        background <- ifelse(is.na(backg),1,backg)
        } else if (identical(backgType,"hybridSD")) {
        background <- ifelse(is.na(backg), pnorm(-2), pnorm(-backg))
        } else if (identical(backgType,"absolute") & 
                 (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
        background <- ifelse(is.na(backg), 
                             pnorm(-2),
                             pnorm((backg-f0)/sqrt(summary(object)$resVar)))
        } else {
        background <- ifelse(is.na(backg),0.1,backg)}
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
                            excess = background - bmr * (1 - background), 
                            additional = background - bmr, 
                            point = bmr)
      } 
      typeVal <- "absolute"
    }
    if (identical(respType, "binomial") & (def %in% c("relative","added","extra", "hybridExc","hybridAdd"))) {
      stop(paste("\"",def, "\" is not available for quantal data", sep=""))
    }
    if (def %in% c("excess","additional") & (backgType %in% c("hybridSD", "hybridPercentile"))) {
      stop(paste("\"",def, "\" does not work with backgType = \"", backgType,"\"", sep=""))
    }
    if (identical(respType, "continuous")) {
      if(identical(slope,"increasing" )) {
        bmrScaled0 <- switch(def, relative = bmr * background + background, 
                             added = bmr + background,
                             point = bmr,
                             extra = bmr*abs(diff(predict(object, data.frame(c(0, Inf)))))
                                     + background,
                             hybridAdd = sqrt(summary(object)$resVar) * 
                                   (qnorm(1 - background) - qnorm(1 - (background + bmr))) + 
                                   f0,
                             hybridExc = sqrt(summary(object)$resVar) * 
                                  (qnorm(1 - background) - qnorm(1 + background - (1 - background)*bmr)) + 
                                   f0)
      } else {
        bmrScaled0 <- switch(def, relative = background - bmr * background, 
                             added = background - bmr,
                             point = bmr,
                             extra = background - bmr*abs(diff(predict(object, data.frame(c(0, Inf))))),
                             hybridAdd = sqrt(summary(object)$resVar) * 
                               (qnorm(background) - qnorm(background + bmr)) + 
                               f0,
                             hybridExc = sqrt(summary(object)$resVar) * 
                               (qnorm(background) - qnorm(background*(1+bmr))) + 
                               f0)
      }
        bmrScaled <- as.numeric(format(bmrScaled0, digits = 5))
        typeVal <- "absolute"
    }
    if (identical(respType, "continuous") & (def %in% c("excess", "additional"))) {
      stop(paste("\"",def, "\" is not available for continuous data", sep=""))
    }
    if (display) {
        cat("Effective response level: ", bmrScaled)
    }
    resMat <- ED(object, bmrScaled, interval = interval, 
            level = 0.9, type = typeVal, display = FALSE)[, 
            c(1, 3), drop = FALSE]
    
    colnames(resMat) <- c("BMD", "BMDL")
    rownames(resMat) <- c("")
    cat("\n\n")
    resMat
    }



