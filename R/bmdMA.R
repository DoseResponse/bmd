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




