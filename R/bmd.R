bmd<-function (object, bmr, backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
               backg=NA, controlSD=NA,
               def = c("excess", "additional", 
                       "relative", "extra", "added", "hybridExc", "hybridAdd", "point"), 
              interval = c("delta", "inv", "profile", "profileGrid"), sandwich.vcov=FALSE, display = TRUE, level=0.95, profileGridSize, profileProgressInfo = TRUE) 
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
  
  # SINGLE CURVE
  if(length(unique(object$dataList$curveid)) == 1){
    slope <- drop(ifelse(object$curve[[1]](0)-object$curve[[1]](Inf)>0,"decreasing","increasing"))
    if(is.na(object$curve[[1]](0)-object$curve[[1]](Inf))){
      slope <- drop(ifelse(object$curve[[1]](0.00000001)-object$curve[[1]](100000000)>0,"decreasing","increasing"))
    }
    if( identical(slope,"increasing" )) {
    f0 <- ifelse(!is.na(coef(object)["c:(Intercept)"]), coef(object)["c:(Intercept)"], object$curve[[1]](0))
    } else {
    f0 <- ifelse(!is.na(coef(object)["d:(Intercept)"]), coef(object)["d:(Intercept)"], object$curve[[1]](0))
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
    
    interval <- match.arg(interval)
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
      resMat <- invBmd(object, bmr, level = 1-2*(1-level), slope=slope, backgType=backgType,
                             backg=backg, catLev=NA, extFactor=10, def=def, useSD=useSD, sandwich.vcov=sandwich.vcov)[, 
                                                                                               c("BMD","BMDL"), drop = FALSE]
      colnames(resMat) <- c("BMD", "BMDL")
      rownames(resMat) <- c("")
      bmdInterval <- invBmd(object, bmr, level = 1-2*(1-level), slope=slope, backgType=backgType,
                                  backg=backg, catLev=NA, extFactor=10, def=def, useSD=useSD, sandwich.vcov=sandwich.vcov)[, 
                                                                                                    c("BMDL", "BMDU"), drop = FALSE]
      colnames(bmdInterval) <- c("Lower CI", "Upper CI")
      rownames(bmdInterval) <- c("")
    }
    
    if(interval == "profile"){
      tmpVals <- ED(object, bmrScaled, interval = "delta", 
                    level = 1-2*(1-level), type = typeVal, vcov. = vcov, display = FALSE)[,c("Estimate", "Lower", "Upper"), drop = FALSE]
      if(backgType %in% c("modelBased", "absolute") & def %in% c("excess", "additional","relative", "extra", "point") 
         & substr(object$fct$name,1,2) %in% c("LL", "LN", "W1", "W2")){    
        if(missing(profileGridSize)){profileGridSize <- 20}
        tmpInterval <- bmdProfileCI(object, bmr = bmr, backgType = backgType, def = def, level = level, gridSize = profileGridSize, 
                                     bmdEst = tmpVals[,"Estimate"], lower = tmpVals[,"Lower"], upper = tmpVals[,"Upper"], slope = slope)
      } else {
        cat("\nReparametrised model not available for chosen backgType, def and dose-response curve. Proceeding with grid search.\n")
        if(missing(profileGridSize)){profileGridSize <- 50}
        tmpInterval <- bmdProfileCIgrid(object, bmr = bmr, backgType = backgType, def = def, level = level, 
                                        gridSize = profileGridSize, progressInfo = profileProgressInfo)
      }
      
      resMat <- matrix(c(tmpVals[,"Estimate"], tmpInterval[1]), nrow = 1)
      colnames(resMat) <- c("BMD", "BMDL")
      rownames(resMat) <- c("")
      
      bmdInterval <- matrix(tmpInterval, nrow = 1)
      colnames(bmdInterval) <- c("Lower CI", "Upper CI")
      rownames(bmdInterval) <- c("")
    }
    
    if(interval == "profileGrid"){
      if(missing(profileGridSize)){profileGridSize <- 50}
      tmpVals <- ED(object, bmrScaled, interval = "delta", 
                    level = 1-2*(1-level), type = typeVal, vcov. = vcov, display = FALSE)[,c("Estimate", "Lower", "Upper"), drop = FALSE]
      
      tmpInterval <- bmdProfileCIgrid(object, bmr = bmr, backgType = backgType, def = def, level = level, 
                                      gridSize = profileGridSize, progressInfo = profileProgressInfo)
      
      resMat <- matrix(c(tmpVals[,"Estimate"], tmpInterval[1]), nrow = 1)
      colnames(resMat) <- c("BMD", "BMDL")
      rownames(resMat) <- c("")
      
      bmdInterval <- matrix(tmpInterval, nrow = 1)
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
  }
  
  # MULTIPLE CURVES FITTED
  if(length(unique(object$dataList$curveid)) > 1){
    
    if(is.null(object$objList)){
      # CURVES ARE NOT FITTED INDEPENDENTLY
      nCurves <- length(unique(object$dataList$curveid))
      
      # Initialise result matrices
      resMat <- matrix(NA, nrow = nCurves, 2)
      bmrScaled <- matrix(NA, nrow = nCurves, 1)
      bmdInterval <- matrix(NA, nrow = nCurves, 2)
      bmdSE <- matrix(NA, nrow = nCurves, 1)
      
      # Compute results for all curves
      def <- match.arg(def)
      interval <- match.arg(interval)
      
      for(iCurve in 1:nCurves){
        slope <- drop(ifelse(object$curve[[1]](0)[,iCurve]-object$curve[[1]](Inf)[,iCurve]>0,"decreasing","increasing"))
        if(is.na(object$curve[[1]](0)[,iCurve]-object$curve[[1]](Inf)[,iCurve])){
          slope <- drop(ifelse(object$curve[[1]](0.00000001)[,iCurve]-object$curve[[1]](100000000)[,iCurve]>0,"decreasing","increasing"))
        }
        
        f0 <- object$curve[[1]](0)[,iCurve]
        
        if(def %in% c("hybridAdd","hybridExc")){
          useSD <- ifelse(!is.na(controlSD),controlSD,sqrt(summary(object)$resVar))
        }
        
        if(!(length(backg) %in% c(1, nCurves))){
          warning('"backg" required to be  of length 1 or number of fitted curves')
          break
        } else if(length(backg) == 1){
          backg <- rep(backg, nCurves)
        }
        if(identical(slope,"increasing" )) {
          if (identical(backgType,"modelBased")) {
            background <- f0
          } else if (identical(backgType,"absolute") & !(def %in% c("relative","extra"))) {
            background <- backg[iCurve]
          } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
            background <- ifelse(is.na(backg[iCurve]),0,backg[iCurve])
          } else {
            background <- ifelse(is.na(backg[iCurve]),1-0.9,1-backg[iCurve])
          }
          if (identical(backgType,"hybridSD")) {
            background <- ifelse(is.na(backg[iCurve]), 1-pnorm(2), 1-pnorm(backg[iCurve]))
          } 
          if (identical(backgType,"absolute") & 
              (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
            background <- ifelse(is.na(backg[iCurve]), 
                                 1-pnorm(2),
                                 1-pnorm((backg[iCurve]-f0)/useSD))
          } 
        } 
        if(identical(slope,"decreasing" )) {
          if (identical(backgType,"modelBased")) {
            background <- f0
          } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
            background <- ifelse(is.na(backg[iCurve]),1,backg[iCurve])
          } else {
            background <- ifelse(is.na(backg[iCurve]),0.1,backg[iCurve])
          }
          if (identical(backgType,"hybridSD")) {
            background <- ifelse(is.na(backg[iCurve]), pnorm(-2), pnorm(-backg[iCurve]))
          } 
          if (identical(backgType,"absolute") & 
              (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
            background <- ifelse(is.na(backg[iCurve]), 
                                 pnorm(-2),
                                 pnorm((backg[iCurve]-f0)/useSD))
          }
        }
        
        respType <- object$type
        if (identical(respType, "binomial")) {
          if(identical(slope,"increasing" )) {
            bmrScaled[iCurve,] <- switch(def, 
                                excess = bmr * (1 - background) + background, 
                                additional = bmr + background, 
                                point = bmr)
          } else {
            bmrScaled[iCurve,] <- switch(def, 
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
            bmrScaled[iCurve,] <- switch(def, 
                                relative = bmr * background + background,
                                added = bmr + background,
                                extra = bmr*abs(diff(predict(object, data.frame(c(0, Inf))))) + background,
                                point = bmr)
          } else {
            bmrScaled[iCurve,] <- switch(def, 
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
                                 extra = bmr*abs(object$curve[[1]](0)[,iCurve]-object$curve[[1]](Inf)[,iCurve])
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
                                 extra = background - bmr*abs(object$curve[[1]](0)[,iCurve]-object$curve[[1]](Inf)[,iCurve]),
                                 hybridAdd = useSD * 
                                   (qnorm(background) - qnorm(background + bmr)) + 
                                   f0,
                                 hybridExc = useSD * 
                                   (qnorm(background) - qnorm(bmr - background*bmr + background)) + 
                                   f0)
          }
          bmrScaled[iCurve,] <- as.numeric(format(bmrScaled0, digits = 5))
          typeVal <- "absolute"
        }
        
        if (identical(respType, "continuous") & (def %in% c("excess", "additional"))) {
          stop(paste("\"",def, "\" is not available for continuous data", sep=""))
        }
        
        if(interval=="delta"){
          if(sandwich.vcov==TRUE){
            resMat[iCurve,] <- ED.bmd(object, bmrScaled[iCurve,], interval = "delta", 
                         level = 1-2*(1-level), type = typeVal, vcov. = sandwich, display = FALSE)[iCurve, 
                                                                                                   c("Estimate", "Lower"), drop = FALSE]
          }
          if(sandwich.vcov==FALSE){
            resMat[iCurve,] <- ED.bmd(object, bmrScaled[iCurve,], interval = interval, 
                         level = 1-2*(1-level), type = typeVal, vcov. = vcov, display = FALSE)[iCurve, 
                                                                                               c("Estimate", "Lower"), drop = FALSE]
          }
          
          
          if(sandwich.vcov==TRUE){
            bmdInterval[iCurve,] <- ED.bmd(object, bmrScaled[iCurve,], interval = interval, 
                              level = 1-2*(1-level), type = typeVal, vcov. = sandwich, display = FALSE)[iCurve, 
                                                                                                        c("Lower", "Upper"), drop = FALSE]
          }
          if(sandwich.vcov==FALSE){
            bmdInterval[iCurve,] <- ED.bmd(object, bmrScaled[iCurve,], interval = interval, 
                              level = 1-2*(1-level), type = typeVal, vcov. = vcov, display = FALSE)[iCurve, 
                                                                                                    c("Lower", "Upper"), drop = FALSE]
          }
        }
        
        if(interval == "inv" | interval == "profile" | interval == "profileGrid"){
          stop(paste("interval \"",inv, "\" is not available for multiple fitted dose-response curves data", sep=""))
        }
        
        if(sandwich.vcov==TRUE){
          bmdSE[iCurve,1] <- ifelse(identical(interval,"delta"),
                               ED.bmd(object, bmrScaled[iCurve,], interval = "delta", 
                                  level = 1-2*(1-level), type = typeVal, vcov. = sandwich, 
                                  display = FALSE)[, c("Std. Error"), drop = FALSE],
                               NA)
        }
        if(sandwich.vcov==FALSE){
          bmdSE[iCurve,1] <- ifelse(identical(interval,"delta"),
                               ED.bmd(object, bmrScaled[iCurve,], interval = "delta", 
                                  level = 1-2*(1-level), type = typeVal, vcov. = vcov, 
                                  display = FALSE)[, c("Std. Error"), drop = FALSE],
                               NA)
        }
      }
      
      rownames(resMat) <- levels(object$dataList$curveid)
      rownames(bmrScaled) <- levels(object$dataList$curveid)
      rownames(bmdInterval) <- levels(object$dataList$curveid)
      rownames(bmdSE) <- levels(object$dataList$curveid)
      
      colnames(resMat) <- c("BMD", "BMDL")
      colnames(bmdSE) <- c("Std. Error")
      colnames(bmdInterval) <- c("Lower CI", "Upper CI")
    } else {
      # CURVES ARE FITTED INDEPENDENTLY
      bmdCall <- function(object){
        bmd(object, bmr, backgType, backg, controlSD,
            def, interval, sandwich.vcov, display = FALSE, level, profileGridSize, profileProgressInfo)
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



