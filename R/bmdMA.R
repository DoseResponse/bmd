bmdMA <- function(modelList, modelWeights, bmr, 
                  backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
                  backg=NA, 
                  def = c("excess", "additional", 
                          "relative", "extra", "added", "hybridExc", "hybridAdd", "point"),
                  respTrans = c("none", "log", "sqrt"),
                  interval = "delta",
                  type = c("curve","bootstrap","Kang","Buckland"),
                  bootstrapType = "nonparametric",
                  R=1000,
                  bootInterval = "percentile",
                  level=0.95,
                  stackingSeed = NULL, stackingSplits = 2,
                  display=TRUE, progressInfo = TRUE){
  nCurves <- ncol(modelList[[1]]$parmMat)
  bmdList<-lapply(modelList, FUN=function(object){bmd(object, bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans,
                                                      interval = interval, display=FALSE, level=level)})  
  if(sum(type == c("curve","bootstrap", "Bootstrap","Kang","Buckland")) != 1){
    cat('Specify model averaging type. Options are "curve", "bootstrap", "Kang" and "Buckland"\n')
  }
  
  if(nCurves == 1){
    # Estimate weights
    if(identical(modelWeights,"AIC")){
      modelWeights0<-exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC))))/
        sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))))
    } else if(identical(modelWeights,"BIC")){
      modelWeights0<-exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC))))/
        sum(exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))))
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
          bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric")
        } else if(identical(bootstrapType,"semiparametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric")
        } else if(identical(bootstrapType,"parametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric")
        }
        
        # if(!oldbootstrap){
          bmdMAboot <- function(data){
            bootModelList <- lapply(modelList, function(model) try(
              eval(substitute(drm(formula = formula0, data = data, fct = model$fct, weights = weights0),
                              list(formula0 = model$call$formula, 
                                   weights0 = model$call$weights))),
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
        # }
        
        # if(oldbootstrap){
        #   drmModelListTmp <-list()
        #   for(i in 1:length(modelList)){
        #     drmModelListTmp[[i]] <- which(!sapply(lapply(bootData, function(x){
        #       try(drm(modelList[[i]]$call$formula, data = x, fct = modelList[[i]][["fct"]]),TRUE)
        #     }
        #     ),function(x) class(x)=="drc"))
        #   }
        #   
        #   non.convergence<-unique(unlist(drmModelListTmp))
        #   if(length(non.convergence)>0){
        #     bootData<-bootData[-non.convergence]
        #   }
        #   bootModelList <-list()
        #   for(i in 1:length(modelList)){
        #     bootModelList[[i]] <- lapply(bootData, function(x){
        #       suppressWarnings(
        #         eval(substitute(drm(formula, data = x, fct = fct0),
        #                         list(formula = modelList[[i]]$call$formula,
        #                              fct0 = modelList[[i]][["fct"]]))
        #         ) # Fitting models using substitute is necessary for Stacking weights
        #       )
        #     }
        #     )
        #   }
        #   
        #   bootModelListTrans <- lapply(1:length(bootModelList[[1]]), function(i) lapply(bootModelList, "[[", i))
        #   
        #   # Compute weights on boot samples
        #   if(identical(modelWeights,"AIC")){
        #     AICList <-suppressWarnings(lapply(bootData, function(x) sapply(modelList, function(y) {AIC(my.fun(x,y))})))
        #     AICtmp <- do.call(rbind,AICList)
        #     modelWeights0 <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
        #   } else if(identical(modelWeights,"BIC")){
        #     BICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
        #     BICtmp <- do.call(rbind,BICList)
        #     modelWeights0 <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
        #   } else if(identical(modelWeights, "Stack")){  
        #     StackList <- lapply(bootModelListTrans, function(x) getStackingWeights(x, stackingSplits))
        #     modelWeights0 <- do.call(cbind, StackList)
        #   } else {
        #     modelWeights0 <- do.call(cbind,rep(list(modelWeights),R))
        #   }
        #   
        #   # Estimate BMD in all models on boot samples
        #   bootbmdList<-list()
        #   for(i in 1:length(modelList)){
        #     bootbmdList[[i]] <- lapply(bootModelList[[i]], function(bootMod){
        #       try(bmd(bootMod,
        #               bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans, interval = interval, 
        #               display=FALSE, level=level)$Results[1,1], silent = TRUE)
        #     }
        #     )
        #   }
        #   bootbmdErrorList <- list()
        #   for(i in 1:length(modelList)){
        #     bootbmdErrorList[[i]] <- which(!sapply(bootbmdList[[i]],function(x) class(x)=="numeric"))
        #   }
        #   
        #   bmd.non.convergence<-unique(unlist(bootbmdErrorList))
        #   if(length(bmd.non.convergence) == length(bootData)){ 
        #     maBMDL <- NA 
        #     maBMDU <- NA
        #   } else {
        #     if(length(bmd.non.convergence)>0){
        #       for(i in 1:length(modelList)){
        #         bootbmdList[[i]]<-bootbmdList[[i]][-bmd.non.convergence]
        #         modelWeights0 <- modelWeights0[,-bmd.non.convergence]
        #       }
        #     }
        #     
        #     boot<-diag(t(matrix(unlist(bootbmdList), ncol = R - length(bmd.non.convergence), byrow = TRUE)) %*% modelWeights0)
        #     boot0<-boot[!is.na(boot)]
        #     
        #     if(bootInterval %in% c("percentile","Percentile")){
        #       maBMDL <- quantile(boot0,p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
        #       maBMDU <- quantile(boot0,p=c(level), na.rm = TRUE)
        #     }
        #     if(identical(bootInterval,"BCa")){
        #       jackData <- list()
        #       for(i in 1:(dim(modelList[[1]]$data)[1])){
        #         jackData[[i]] <- modelList[[1]]$data[-i,]
        #       }
        #       bootJackList <-list()
        #       for(i in 1:length(modelList)){
        #         bootJackList[[i]] <- sapply(jackData, function(x){
        #           suppressWarnings(bmd(drm(modelList[[i]]$call$formula, data = x, type = modelList[[i]]$type, 
        #                                    fct = modelList[[i]][["fct"]]),
        #                                bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans, interval = interval, 
        #                                display=FALSE, level=level)$Results[1])
        #         }
        #         )
        #       }
        #       if(identical(modelWeights,"AIC")){
        #         AICJackList <-list()
        #         for(i in 1:length(modelList)){
        #           AICJackList[[i]] <- sapply(jackData, function(x){
        #             suppressWarnings(AIC(drm(modelList[[i]]$call$formula, data = x, type = modelList[[i]]$type, fct = modelList[[i]][["fct"]])))
        #           }
        #           )
        #         }
        #         AICtmp <- do.call(rbind,AICJackList)
        #         modelWeightsJack <- t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))
        #       } else if(identical(modelWeights,"BIC")){
        #         BICJackList <-list()
        #         for(i in 1:length(modelList)){
        #           BICJackList[[i]] <- sapply(jackData, function(x){
        #             suppressWarnings(BIC(drm(modelList[[i]]$call$formula, data = x, type = modelList[[i]]$type, fct = modelList[[i]][["fct"]])))
        #           }
        #           )
        #         }
        #         BICtmp <- do.call(rbind,BICJackList)
        #         modelWeightsJack <- t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))
        #       } else {
        #         modelWeightsJack <- do.call(cbind,rep(list(modelWeights),dim(modelList[[1]]$data)[1]))
        #       }
        #       
        #       bootjack<-diag(t(do.call(rbind,bootJackList)) %*% modelWeightsJack)
        #       
        #       maBMDL <- as.numeric(BCa(obs = maBMD, data = modelList[[1]]$data, boot0, bootjack, level = level)[1])
        #       maBMDU <- "Not available for BCa bootstrap"
        #     }
        #   }
        # }
      }
      
      if(identical(type,"curve")){
        bmrScaled0<-sum(sapply(bmdList, function(x){x$bmrScaled})*modelWeights0)
        maBMD <- bmdMACurve(modelList,modelWeights0,bmrScaled0)$Results[1]
        
        # Bootstrap
        if(identical(bootstrapType,"nonparametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric")
        } else if(identical(bootstrapType,"semiparametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric")
        } else if(identical(bootstrapType,"parametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric")
        }
        
        # if(!oldbootstrap){
          bmdMACurveboot <- function(data){
            bootModelList <- lapply(modelList, function(model){try(
              eval(substitute(
                drm(formula = formula0, data = data, fct = model$fct, weights = weights0, type = model$type),
                list(formula0 = model$call$formula, 
                     weights0 = model$call$weights))),
              silent = TRUE)
            })
            
            modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
            
            bootModelList <- bootModelList[!modelConvergenceError]
            bootBmdList <- lapply(bootModelList, 
                                  function(object){
                                    bmd(object, bmr = bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans,
                                        interval = "delta", display=FALSE)})
            
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
            bootBmrScaled0 <- sum(sapply(bootBmdList, function(x){x$bmrScaled})*bootModelWeights0)
            
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
        # }
        
        # if(oldbootstrap){
        #   drmModelListTmp <-list()
        #   for(i in 1:length(modelList)){
        #     drmModelListTmp[[i]] <- which(!sapply(lapply(bootData, function(x){
        #       try(drm(modelList[[i]]$call$formula, data = x, fct = modelList[[i]][["fct"]]),TRUE)
        #     }
        #     ),function(x) class(x)=="drc"))
        #   }
        #   
        #   non.convergence<-unique(unlist(drmModelListTmp))
        #   if(length(non.convergence)>0){
        #     bootData<-bootData[-non.convergence]
        #   }
        #   bootModelList <-list()
        #   for(i in 1:length(modelList)){
        #     bootModelList[[i]] <- lapply(bootData, function(x){
        #       suppressWarnings(
        #         eval(substitute(drm(formula, data = x, fct = fct0),
        #                         list(formula = modelList[[i]]$call$formula,
        #                              fct0 = modelList[[i]][["fct"]]))
        #         ) # Fitting models using substitute is necessary for Stacking weights
        #       )
        #     }
        #     )
        #   }
        #   
        #   bootModelListTrans <- lapply(1:length(bootModelList[[1]]), function(i) lapply(bootModelList, "[[", i))
        #   
        #   if(identical(modelWeights,"AIC")){
        #     AICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(AIC(my.fun(x,y)))))
        #     AICtmp <- do.call(rbind,AICList)
        #     modelWeights0 <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
        #   } else if(identical(modelWeights,"BIC")){
        #     BICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
        #     BICtmp <- do.call(rbind,BICList)
        #     modelWeights0 <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
        #   } else if(identical(modelWeights, "Stack")){  
        #     modelWeightsList <- lapply(bootModelListTrans, function(x) getStackingWeights(x, stackingSplits))
        #   } else {
        #     modelWeights0 <- do.call(cbind,rep(list(modelWeights),R))
        #   }
        #   
        #   if(!identical(modelWeights, "Stack")){ 
        #     modelWeightsList <- lapply(1:ncol(modelWeights0),function(i) modelWeights0[,i])
        #   }
        #   
        #   bootbmrList<-list()
        #   for(i in 1:length(modelList)){
        #     bootbmrList[[i]] <- sapply(bootModelList[[i]], function(bootMod){
        #       as.numeric(try(bmd(bootMod,
        #                          bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans, interval = interval, 
        #                          display=FALSE, level=level)$bmrScaled, silent = TRUE))
        #     }
        #     )
        #   }
        #   bootbmrListTrans <- lapply(1:length(bootbmrList[[1]]), function(i) sapply(bootbmrList, "[[", i))
        #   
        #   LLimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]]))[2]/10000
        #   ULimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]],decreasing=TRUE))[1]
        #   funk<-function(x,y,z){try(bmdMACurve(x,y,z,searchInterval=c(LLimit,ULimit))$Results[1],TRUE)}
        #   bmrScaledList<-as.list(rowSums(do.call(rbind,modelWeightsList)*do.call(rbind,bootbmrListTrans)))
        #   
        #   boot<-mapply(funk,bootModelListTrans,modelWeightsList,bmrScaledList)
        #   boot0<-suppressWarnings(as.numeric(boot[!is.na(as.numeric(boot))]))
        #   
        #   if(bootInterval %in% c("percentile","Percentile")){
        #     maBMDL <- quantile(boot0,p=c(1-level), na.rm = FALSE) # ABC percentile lims.  
        #     maBMDU <- quantile(boot0,p=c(level), na.rm = FALSE)
        #   }
        #   if(identical(bootInterval,"BCa")){
        #     jackData <- list()
        #     for(i in 1:(dim(modelList[[1]]$data)[1])){
        #       jackData[[i]] <- modelList[[1]]$data[-i,]
        #     }
        #     
        #     bootJackModelList <- lapply(jackData, function(x) lapply(modelList, function(y) suppressWarnings(my.fun(x,y))))
        #     
        #     
        #     if(identical(modelWeights,"AIC")){
        #       AICJackList <-lapply(bootJackModelList, function(x) sapply(x, AIC))#lapply(jackData, function(x) sapply(modelList, function(y) suppressWarnings(AIC(my.fun(x,y)))))
        #       AICtmp <- do.call(rbind,AICJackList)
        #       modelWeightsJack <- t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))
        #     } else if(identical(modelWeights,"BIC")){
        #       BICJackList <-lapply(bootJackModelList, function(x) sapply(x, BIC))#lapply(jackData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
        #       BICtmp <- do.call(rbind,BICJackList)
        #       modelWeightsJack <- t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))
        #     } else {
        #       modelWeightsJack <- do.call(cbind,rep(list(modelWeights),dim(modelList[[1]]$data)[1]))
        #     }
        #     modelWeightsJackList <- lapply(1:ncol(modelWeightsJack),function(i) modelWeightsJack[,i])
        #     modelWeightsJackListTrans <- lapply(1:length(modelWeightsJackList[[1]]), function(i) sapply(modelWeightsJackList, "[[", i))
        #     
        #     jackbmrList<-list()
        #     for(i in 1:length(modelList)){
        #       jackbmrList[[i]] <- sapply(jackData, function(x){
        #         suppressWarnings(bmd(drm(modelList[[i]]$call$formula, data = x, type = modelList[[i]]$type, 
        #                                  fct = modelList[[i]][["fct"]]),
        #                              bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans, interval = interval, 
        #                              display=FALSE, level=level)$bmrScaled)
        #       }
        #       )
        #     }
        #     
        #     jackbmrListTrans <- lapply(1:length(jackbmrList[[1]]), function(i) sapply(jackbmrList, "[[", i))
        #     
        #     LLimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]]))[2]/10000
        #     ULimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]],decreasing=TRUE))[1]
        #     funk<-function(x,y,z){
        #       as.numeric(try(bmdMACurve(x,y,z,searchInterval=c(LLimit,ULimit))$Results[1], silent = TRUE))
        #     }
        #     
        #     bmrScaledJackTrans<-as.list(diag(do.call(cbind,modelWeightsJackList) %*% do.call(cbind,jackbmrListTrans)))
        #     
        #     bootjack<-mapply(funk,bootJackModelList,modelWeightsJackListTrans,bmrScaledJackTrans)
        #     
        #     maBMDL <- as.numeric(BCa(obs = maBMD, data = modelList[[1]]$data, boot0, bootjack, level = level)[1])
        #     maBMDU <- "Not available for BCa bootstrap"
        #   }
        # }
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
          bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric",aggregated = FALSE)
        } else if(identical(bootstrapType,"semiparametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric",aggregated = FALSE)
        } else if(identical(bootstrapType,"parametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric",aggregated = FALSE)
        }
        
        # if(!oldbootstrap){
          bmdMAboot <- function(data){
            bootModelList <- lapply(modelList, function(model) try(
              eval(substitute(drm(formula = formula0, data = data, fct = model$fct, weights = weights0, type = "binomial"),
                              list(formula0 = model$call$formula, 
                                   weights0 = model$call$weights))),
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
              
              maBMDL <- as.numeric(BCa(obs = maBMD, data = modelList[[1]]$data, boot0, bootjack, level = level)[1])
              maBMDU <- "Not available for BCa bootstrap"
            }
          }
        # }
        
        # if(oldbootstrap){
        #   drmModelListTmp <-list()
        #   for(i in 1:length(modelList)){
        #     drmModelListTmp[[i]] <- which(!sapply(lapply(bootData, function(x){
        #       try(drm(modelList[[i]]$call$formula, data = x, type="binomial", fct = modelList[[i]][["fct"]]),TRUE)
        #     }
        #     ),function(x) class(x)=="drc"))
        #   }
        #   
        #   non.convergence<-unique(unlist(drmModelListTmp))
        #   if(length(non.convergence)>0){
        #     bootData<-bootData[-non.convergence]
        #   }
        #   bootModelList <-list()
        #   for(i in 1:length(modelList)){
        #     bootModelList[[i]] <- lapply(bootData, function(x){
        #       suppressWarnings(
        #         eval(substitute(drm(formula, data = x, fct = fct0, type = "binomial"),
        #                         list(formula = modelList[[i]]$call$formula,
        #                              fct0 = modelList[[i]][["fct"]]))
        #         ) # Fitting models using substitute is necessary for Stacking weights
        #       )
        #     }
        #     )
        #   }
        #   
        #   bootModelListTrans <- lapply(1:length(bootModelList[[1]]), function(i) lapply(bootModelList, "[[", i))
        #   
        #   # Compute weights
        #   if(identical(modelWeights,"AIC")){
        #     AICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(AIC(my.fun(x,y)))))
        #     AICtmp <- do.call(rbind,AICList)
        #     modelWeights0 <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
        #   } else if(identical(modelWeights,"BIC")){
        #     BICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
        #     BICtmp <- do.call(rbind,BICList)
        #     modelWeights0 <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
        #   } else if(identical(modelWeights, "Stack")){  
        #     StackList <- lapply(bootModelListTrans, function(x) getStackingWeights(x, stackingSplits))
        #     modelWeights0 <- do.call(cbind, StackList)
        #   } else {
        #     modelWeights0 <- do.call(cbind,rep(list(modelWeights),R))
        #   }
        #   
        #   # Estimate BMD in all models on boot samples
        #   bootbmdList<-list()
        #   for(i in 1:length(modelList)){
        #     bootbmdList[[i]] <- lapply(bootModelList[[i]], function(bootMod){
        #       try(bmd(bootMod,
        #               bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans, interval = interval, 
        #               display=FALSE, level=level)$Results[1,1], silent = TRUE)
        #     }
        #     )
        #   }
        #   bootbmdErrorList <- list()
        #   for(i in 1:length(modelList)){
        #     bootbmdErrorList[[i]] <- which(!sapply(bootbmdList[[i]],function(x) class(x)=="numeric"))
        #   }
        #   
        #   bmd.non.convergence<-unique(unlist(bootbmdErrorList))
        #   if(length(bmd.non.convergence) == length(bootData)){ 
        #     maBMDL <- NA 
        #     maBMDU <- NA
        #   } else {
        #     if(length(bmd.non.convergence)>0){
        #       for(i in 1:length(modelList)){
        #         bootbmdList[[i]]<-bootbmdList[[i]][-bmd.non.convergence]
        #       }
        #       modelWeights0 <- modelWeights0[,-bmd.non.convergence]
        #     }
        #     
        #     boot<-diag(t(matrix(unlist(bootbmdList), ncol = length(bootbmdList[[1]]), byrow = TRUE)) %*% modelWeights0)
        #     boot0<-boot[!is.na(boot)]
        #     
        #     if(bootInterval %in% c("percentile","Percentile")){
        #       maBMDL <- quantile(boot0,p=c(1-level), na.rm = FALSE) # ABC percentile lims.  
        #       maBMDU <- quantile(boot0,p=c(level), na.rm = FALSE)
        #     }
        #     if(identical(bootInterval,"BCa")){
        #       data.str <- modelList[[1]]$data
        #       data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
        #       data.e<-expandBinomial(data.str, 
        #                              number = "number",
        #                              total = "weights",
        #                              dose = as.character(modelList[[1]]$call$formula[[3]]))
        #       jackData <- list()
        #       for(i in 1:(dim(data.e)[1])){
        #         jackData[[i]] <- data.e[-i,]
        #       }
        #       bootJackList <-list()
        #       for(i in 1:length(modelList)){
        #         bootJackList[[i]] <- sapply(jackData, function(x){
        #           as.numeric(try(bmd(my.fun2(x,modelList[[i]]),
        #                              bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans, interval = interval, 
        #                              display=FALSE, level=level)$Results[1], silent = TRUE))
        #         }
        #         )
        #       }
        #       if(identical(modelWeights,"AIC")){
        #         AICJackList <-list()
        #         for(i in 1:length(modelList)){
        #           AICJackList[[i]] <- sapply(jackData, function(x){
        #             as.numeric(try(AIC(my.fun2(x,modelList[[i]])), silent = TRUE))
        #           }
        #           )
        #         }
        #         AICtmp <- t(do.call(rbind,AICJackList))
        #         modelWeightsJack <- t(t(exp(-t(AICtmp - do.call(pmin, (as.data.frame(AICtmp))))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
        #       } else if(identical(modelWeights,"BIC")){
        #         BICJackList <-list()
        #         for(i in 1:length(modelList)){
        #           BICJackList[[i]] <- sapply(jackData, function(x){
        #             suppressWarnings(BIC(my.fun2(x,modelList[[i]])))
        #           }
        #           )
        #         }
        #         BICtmp <- t(do.call(rbind,BICJackList))
        #         modelWeightsJack <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
        #       } else {
        #         modelWeightsJack <- do.call(cbind,rep(list(modelWeights),dim(modelList[[1]]$data)[1]))
        #       }
        #       
        #       bootjack<-diag(t(do.call(rbind,bootJackList)) %*% modelWeightsJack)
        #       
        #       maBMDL <- as.numeric(BCa(obs = maBMD, data = data.e, boot0, bootjack, level = level)[1])
        #       maBMDU <- "Not available for BCa bootstrap"
        #     }
        #   }
        # }
      }
      
      if(identical(type,"curve")){
        bmrScaled0<-sum(sapply(bmdList, function(x){x$bmrScaled})*modelWeights0)
        maBMD <- bmdMACurve(modelList,modelWeights0,bmrScaled0)$Results[1]
        
        # Bootstrap
        if(identical(bootstrapType,"nonparametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric",aggregated = FALSE)
        } else if(identical(bootstrapType,"semiparametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric",aggregated=FALSE)
        } else if(identical(bootstrapType,"parametric")){
          bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric",aggregated=FALSE)
        }
        
        # if(!oldbootstrap){
          bmdMACurveboot <- function(data){
            bootModelList <- lapply(modelList, function(model){try(
              eval(substitute(
                drm(formula = formula0, data = data, fct = model$fct, weights = weights0, type = model$type),
                list(formula0 = model$call$formula, 
                     weights0 = model$call$weights))),
              silent = TRUE)
            })
            
            modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
            
            bootModelList <- bootModelList[!modelConvergenceError]
            bootBmdList <- lapply(bootModelList, 
                                  function(object){
                                    bmd(object, bmr = bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans,
                                        interval = "delta", display=FALSE)})
            
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
            
            bootBmrScaled0 <- sum(sapply(bootBmdList, function(x){x$bmrScaled})*bootModelWeights0)
            
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
              
              maBMDL <- as.numeric(BCa(obs = maBMD, data = modelList[[1]]$data, boot0, bootjack, level = level)[1])
              maBMDU <- "Not available for BCa bootstrap"
            }
          }
        # }
        
        # if(oldbootstrap){
        #   drmModelListTmp <-list()
        #   for(i in 1:length(modelList)){
        #     drmModelListTmp[[i]] <- which(!sapply(lapply(bootData, function(x){
        #       try(drm(modelList[[i]]$call$formula, data = x, type="binomial", fct = modelList[[i]][["fct"]]),TRUE)
        #     }
        #     ),function(x) class(x)=="drc"))
        #   }
        #   
        #   non.convergence<-unique(unlist(drmModelListTmp))
        #   if(length(non.convergence)>0){
        #     bootData<-bootData[-non.convergence]
        #   }
        #   
        #   bootModelList <-list()
        #   for(i in 1:length(modelList)){
        #     bootModelList[[i]] <- lapply(bootData, function(x){
        #       suppressWarnings(
        #         eval(substitute(drm(formula, data = x, fct = fct0, type = "binomial"),
        #                         list(formula = modelList[[i]]$call$formula,
        #                              fct0 = modelList[[i]][["fct"]]))
        #         )
        #       )
        #     }
        #     )
        #   }
        #   
        #   bootModelListTrans <- lapply(1:length(bootModelList[[1]]), function(i) lapply(bootModelList, "[[", i))
        #   
        #   if(identical(modelWeights,"AIC")){
        #     AICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(AIC(my.fun(x,y)))))
        #     AICtmp <- do.call(rbind,AICList)
        #     modelWeights0 <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
        #     
        #   } else if(identical(modelWeights,"BIC")){
        #     BICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
        #     BICtmp <- do.call(rbind,BICList)
        #     modelWeights0 <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
        #   } else if(identical(modelWeights, "Stack")){  
        #     modelWeightsList <- lapply(bootModelListTrans, function(x) getStackingWeights(x, stackingSplits))
        #   } else {
        #     modelWeights0 <- do.call(cbind,rep(list(modelWeights),R))
        #   }
        #   
        #   if(!identical(modelWeights, "Stack")){ 
        #     modelWeightsList <- lapply(1:ncol(modelWeights0),function(i) modelWeights0[,i])
        #   }
        #   
        #   bootbmrList<-list()
        #   for(i in 1:length(modelList)){
        #     bootbmrList[[i]] <- lapply(bootModelList[[i]], function(bootMod){
        #       try(bmd(bootMod,
        #               bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans, interval = interval, 
        #               display=FALSE, level=level)$bmrScaled, silent = TRUE)
        #     }
        #     )
        #   }
        #   bootbmrListTrans <- lapply(1:length(bootbmrList[[1]]), function(i) sapply(bootbmrList, "[[", i))
        #   
        #   bootbmrErrorList <- list()
        #   for(i in 1:length(modelList)){
        #     bootbmrErrorList[[i]] <- which(sapply(bootbmrList[[i]],function(x) is.na(as.numeric(x)))) #class(x)=="numeric"
        #   }
        #   
        #   bmr.non.convergence<-unique(unlist(bootbmrErrorList))
        #   if(length(bmr.non.convergence) > 0){
        #     bootbmrListTrans <- bootbmrListTrans[-bmr.non.convergence]
        #     modelWeightsList <- modelWeightsList[-bmr.non.convergence]
        #     bootModelListTrans <- bootModelListTrans[-bmr.non.convergence]
        #   }
        #   
        #   LLimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]]))[2]/10000
        #   ULimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]],decreasing=TRUE))[1]
        #   funk<-function(x,y,z){
        #     try(bmdMACurve(x,y,z,searchInterval=c(LLimit,ULimit))$Results[1], silent = TRUE)}
        #   bmrScaledList<-as.list(rowSums(do.call(rbind,modelWeightsList)*do.call(rbind,bootbmrListTrans)))
        #   
        #   boot<-mapply(funk,bootModelListTrans,modelWeightsList,bmrScaledList)
        #   boot0<-suppressWarnings(as.numeric(boot[!is.na(as.numeric(boot))]))
        #   
        #   if(bootInterval %in% c("percentile","Percentile")){
        #     maBMDL <- quantile(boot0,p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
        #     maBMDU <- quantile(boot0,p=c(level), na.rm = TRUE)
        #   }
        #   if(identical(bootInterval,"BCa")){
        #     data.str <- modelList[[1]]$data
        #     data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
        #     data.e<-expandBinomial(data.str, 
        #                            number = "number",
        #                            total = "weights",
        #                            dose = as.character(modelList[[1]]$call$formula[[3]]))
        #     
        #     jackData <- list()
        #     for(i in 1:(dim(data.e)[1])){
        #       jackData[[i]] <- data.e[-i,]
        #     }
        #     
        #     bootJackModelList <- lapply(jackData, function(x) lapply(modelList, function(y) suppressWarnings(my.fun2(x,y))))
        #     
        #     if(identical(modelWeights,"AIC")){
        #       AICJackList <-list()
        #       for(i in 1:length(modelList)){
        #         AICJackList[[i]] <- sapply(jackData, function(x){
        #           as.numeric(try(AIC(my.fun2(x,modelList[[i]])), silent = TRUE))
        #         }
        #         )
        #       }
        #       AICtmp <- t(do.call(rbind,AICJackList))
        #       modelWeightsJack <- t(t(exp(-t(AICtmp - do.call(pmin, (as.data.frame(AICtmp))))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
        #     } else if(identical(modelWeights,"BIC")){
        #       BICJackList <-list()
        #       for(i in 1:length(modelList)){
        #         BICJackList[[i]] <- sapply(jackData, function(x){
        #           suppressWarnings(BIC(my.fun2(x,modelList[[i]])))
        #         }
        #         )
        #       }
        #       BICtmp <- t(do.call(rbind,BICJackList))
        #       modelWeightsJack <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
        #     } else {
        #       modelWeightsJack <- do.call(cbind,rep(list(modelWeights),dim(modelList[[1]]$data)[1]))
        #     }
        #     
        #     # modelWeightsJackList <- lapply(1:ncol(modelWeightsJack),function(i) modelWeightsJack[,i])
        #     # modelWeightsJackListTrans <- lapply(1:length(modelWeightsJackList[[1]]), function(i) sapply(modelWeightsJackList, "[[", i))
        #     
        #     jackbmrList<-list()
        #     for(i in 1:length(modelList)){
        #       jackbmrList[[i]] <- sapply(jackData, function(x){
        #         as.numeric(try(bmd(my.fun2(x,modelList[[i]]),
        #                            bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans, interval = interval, 
        #                            display=FALSE, level=level)$bmrScaled,silent=TRUE))
        #       }
        #       )
        #     }
        #     
        #     jackbmrListTrans <- lapply(1:length(jackbmrList[[1]]), function(i) sapply(jackbmrList, "[[", i))
        #     
        #     LLimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]]))[2]/10000
        #     ULimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]],decreasing=TRUE))[1]
        #     funk<-function(x,y,z){
        #       as.numeric(try(bmdMACurve(x,y,z,searchInterval=c(LLimit,ULimit))$Results[1], silent = TRUE))
        #     }
        #     bmrScaledJack<-as.list(diag(do.call(cbind,modelWeightsJackList) %*% do.call(cbind,jackbmrListTrans)))
        #     
        #     bootjack<-mapply(funk,bootJackModelList,modelWeightsJackListTrans,bmrScaledJack)
        #     
        #     maBMDL <- as.numeric(BCa(obs = maBMD, data = data.e, boot0, bootjack, level = level)[1])
        #     maBMDU <- "Not available for BCa bootstrap"
        #   } 
        # }
      }
    }
  } 
  
  if (nCurves > 1){
    if(is.null(modelList[[1]]$objList)){
      # Estimate weights
      if(identical(modelWeights,"AIC")){
        modelWeights0<-exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC))))/
          sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))))
      } else if(identical(modelWeights,"BIC")){
        modelWeights0<-exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC))))/
          sum(exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))))
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
            bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric")
          } else if(identical(bootstrapType,"semiparametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric")
          } else if(identical(bootstrapType,"parametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric")
          }
          
          bmdMAboot <- function(data){
            data[[as.character(modelList[[1]]$call$curveid)]] <- data[[paste0("orig.", as.character(modelList[[1]]$call$curveid))]]
            if(is.null(modelList[[1]]$call$pmodels)){
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula = formula0, data = data, fct = model$fct, weights = weights0, 
                      curveid = model$call$curveid, type = model$type, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula, 
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid))),
                silent = TRUE)
              })
            } else {
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula0, weights = weights0, curveid = curveid0, pmodels = pmodels0,
                      data = data, type = model$type, fct = model$fct, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula,
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       pmodels0 = model$call$pmodels)
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
                jackData[[i]] <- modelList[[1]]$data[-i,]
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
            bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric")
          } else if(identical(bootstrapType,"semiparametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric")
          } else if(identical(bootstrapType,"parametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric")
          }
          
          bmdMACurveboot <- function(data){
            data[[as.character(modelList[[1]]$call$curveid)]] <- data[[paste0("orig.", as.character(modelList[[1]]$call$curveid))]]
            if(is.null(modelList[[1]]$call$pmodels)){
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula = formula0, data = data, fct = model$fct, weights = weights0, 
                      curveid = curveid0, type = model$type, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula, 
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid))),
                silent = TRUE)
              })
            } else {
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula0, weights = weights0, curveid = curveid0, pmodels = pmodels0,
                      data = data, type = model$type, fct = model$fct, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula,
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       pmodels0 = model$call$pmodels)
                )),
                silent = TRUE)
              })
            }
            modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
            
            bootModelList <- bootModelList[!modelConvergenceError]
            bootBmdList <- lapply(bootModelList, 
                                  function(object){
                                    bmd(object, bmr = bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans,
                                        interval = "delta", display=FALSE)})
            
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
            
            bootBmrScaled0 <- colSums(bootModelWeights0 * t(sapply(bootBmdList, function(x){x$bmrScaled})))
            
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
                jackData[[i]] <- modelList[[1]]$data[-i,]
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
            bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric",aggregated = FALSE)
          } else if(identical(bootstrapType,"semiparametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric",aggregated=FALSE)
          } else if(identical(bootstrapType,"parametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric",aggregated=FALSE)
          }
          
          bmdMAboot <- function(data){
            if(is.null(modelList[[1]]$call$pmodels)){
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula = formula0, data = data, fct = model$fct, weights = weights0, 
                      curveid = curveid0, type = model$type, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula, 
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid))),
                silent = TRUE)
              })
            } else {
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula0, weights = weights0, curveid = curveid0, pmodels = pmodels0,
                      data = data, type = model$type, fct = model$fct, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula,
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       pmodels0 = model$call$pmodels)
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
            bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric",aggregated = FALSE)
          } else if(identical(bootstrapType,"semiparametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric",aggregated=FALSE)
          } else if(identical(bootstrapType,"parametric")){
            bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric",aggregated=FALSE)
          }
          
          bmdMACurveboot <- function(data){
            if(is.null(modelList[[1]]$call$pmodels)){
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula = formula0, data = data, fct = model$fct, weights = weights0, 
                      curveid = curveid0, type = model$type, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula, 
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid))),
                silent = TRUE)
              })
            } else {
              bootModelList <- lapply(modelList, function(model){try(
                eval(substitute(
                  drm(formula0, weights = weights0, curveid = curveid0, pmodels = pmodels0,
                      data = data, type = model$type, fct = model$fct, control = drmc(noMessage = TRUE)),
                  list(formula0 = model$call$formula,
                       weights0 = model$call$weights,
                       curveid0 = model$call$curveid,
                       pmodels0 = model$call$pmodels)
                )),
                silent = TRUE)
              })
            }
            modelConvergenceError <- sapply(bootModelList, function(mod_try) inherits(mod_try, "try-error"))
            
            bootModelList <- bootModelList[!modelConvergenceError]
            bootBmdList <- lapply(bootModelList, 
                                  function(object){
                                    bmd(object, bmr = bmr, backgType = backgType, backg = backg, def = def, respTrans = respTrans,
                                        interval = "delta", display=FALSE)})
            
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
            
            bootBmrScaled0 <- colSums(bootModelWeights0 * t(sapply(bootBmdList, function(x){x$bmrScaled})))
            
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
              
              maBMDL <- sapply(1:nCurves, function(i) as.numeric(BCa(obs = maBMD[i], data = modelList[[1]]$data, boot0[,i], bootjack[,i], level = level)[1]))
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
      maBMDL <- sapply(bmdMAList, function(x) x$Interval[,1])
      maBMDU <- sapply(bmdMAList, function(x) x$Interval[,2])
      maBMDse <- sapply(bmdMAList, function(x) x$SE[,1])
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
               Interval = intMat,
               SE = seMat,
               modelWeights = modelWeights0)
  if(display){
    print(resMat)
  }
  
  class(resBMD) <- "bmd"
  invisible(resBMD)
}




