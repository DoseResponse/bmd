bmd<-function (object, bmr, backgType = c("modelBased","absolute","hybridSD","hybridPercentile"),
               backg=NA, 
               def = c("excess", "additional", 
                       "relative", "extra","absolute", "hybridExc", "hybridAdd", "point"), 
              interval = "delta", display = FALSE) 
{
    slope <- ifelse(as.numeric(predict(object,data.frame(0))-predict(object,data.frame(Inf)))>0,"decreasing","increasing")
    if( identical(slope,"increasing" )) {
    f0 <- ifelse(!is.na(coef(object)["c:(Intercept)"]), coef(object)["c:(Intercept)"], predict(object,data.frame(0)))
    } else {
    f0 <- ifelse(!is.na(coef(object)["d:(Intercept)"]), coef(object)["d:(Intercept)"], predict(object,data.frame(0)))
    }
    if(identical(slope,"increasing" )) {
      if (identical(backgType,"modelBased")) {
        background <- f0
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
        background <- ifelse(is.na(backg),0.9,backg)}
    } else {
      if (identical(backgType,"modelBased")) {
        background <- f0
        } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
        background <- ifelse(is.na(backg),1,backg)
        } else if (identical(backgType,"hybridSD")) {
        background <- ifelse(is.na(backg), 1-pnorm(2), 1-pnorm(backg))
        } else if (identical(backgType,"absolute") & 
                 (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
        background <- ifelse(is.na(backg), 
                             1-pnorm(2),
                             1-pnorm((f0-backg)/sqrt(summary(object)$resVar)))
        } else {
        background <- ifelse(is.na(backg),0.9,backg)}
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
    if (identical(respType, "binomial") & (def %in% c("relative","absolute","hybridExc","hybridAdd"))) {
      stop(paste("\"",def, "\" is not available for quantal data", sep=""))
    }
    if (identical(respType, "continuous")) {
      if(identical(slope,"increasing" )) {
        bmrScaled0 <- switch(def, relative = bmr * background + background, 
                             absolute = bmr + background,
                             point = bmr,
                             std = bmr*sqrt(summary(object)$resVar) + background,
                             extra = bmr*abs(diff(predict(object, data.frame(c(0, Inf)))))
                                     + background,
                             hybridAdd = sqrt(summary(object)$resVar) * 
                                   (qnorm(1 - background) - qnorm(1 - (background + bmr))) + 
                                   coef(object)["c:(Intercept)"],
                             hybridExc = sqrt(summary(object)$resVar) * 
                                  (qnorm(1 - background) - qnorm(1 + background - (1 - background)*bmr)) + 
                                  coef(object)["c:(Intercept)"])
      } else {
        bmrScaled0 <- switch(def, relative = background - bmr * background, 
                             absolute = background - bmr,
                             point = bmr,
                             std = background - bmr*sqrt(summary(object)$resVar),
                             extra = background - bmr*abs(diff(predict(object, data.frame(c(0, Inf))))),
                             hybridAdd = sqrt(summary(object)$resVar) * 
                               (qnorm(1 - background) - qnorm(1 - (background + bmr))) + 
                               coef(object)["c:(Intercept)"],
                             hybridExc = sqrt(summary(object)$resVar) * 
                               (qnorm(1 - background) - qnorm(1 + background - (1 - background)*bmr)) + 
                               coef(object)["c:(Intercept)"])
      }
        bmrScaled <- as.numeric(format(bmrScaled0, digits = 5))
        typeVal <- "absolute"
    }
    if (identical(respType, "continuous") & (def %in% c("excess", "additional"))) {
      stop(paste("\"",def, "\" is not available for continuous data", sep=""))
    }
    if (display) {
        cat("Effective response level: ", bmrScaled, 
            "\n\n")
    }
    resMat <- ED(object, bmrScaled, interval = interval, 
            level = 0.9, type = typeVal, display = FALSE)[, 
            c(1, 3), drop = FALSE]
    
    colnames(resMat) <- c("BMD", "BMDL")
    rownames(resMat) <- c("")
    cat("\n\n")
    resMat
    }
