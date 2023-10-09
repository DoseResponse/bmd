bmdMA <- function(modelList, modelWeights, bmr, 
                  backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
                  backg=NA, 
                  def = c("excess", "additional", 
                          "relative", "extra", "added", "hybridExc", "hybridAdd", "point"),
                  interval = "delta",
                  type = c("curve","bootstrap","Kang","Buckland"),
                  bootstrapType = "nonparametric",
                  R=1000,
                  bootInterval = "percentile",
                  level=0.95,
                  display=TRUE){
  
  bmdList<-lapply(modelList, FUN=function(object){bmd(object, bmr, backgType = backgType, backg = backg, def = def, 
                                                      interval = interval, display=FALSE, level=level)})  
  
  if(identical(modelList[[1]]$type,"continuous")){
    my.fun<-function(x,y){drm(y$call$formula, data = x, fct = y[["fct"]])}
    
    if(identical(modelWeights,"AIC")){
      modelWeights0<-exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC))))/
        sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))))
      exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC))))/
        sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))))
    } else if(identical(modelWeights,"BIC")){
      modelWeights0<-exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC))))/
        sum(exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))))
    } else if(identical(modelWeights,"Stack")){
      modelWeights0<-getStackingWeights(modelList)
    } else {
      modelWeights0 <- modelWeights
    }
    
    if(identical(type,"Kang")){
      maBMD <- sum(modelWeights0 * sapply(bmdList, function(x){x$Results[,1]}))
      maBMDL <- sum(modelWeights0 * sapply(bmdList, function(x){x$Results[,2]}))
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
      if(identical(bootstrapType,"nonparametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric")
      } else if(identical(bootstrapType,"semiparametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric")
      } else if(identical(bootstrapType,"parametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric")
      }
      maBMD <- sum(modelWeights0 * sapply(bmdList, function(x){x$Results[,1]}))
      
      drmModelListTmp <-list()
      for(i in 1:length(modelList)){
        drmModelListTmp[[i]] <- which(!sapply(lapply(bootData, function(x){
          try(drm(modelList[[i]]$call$formula, data = x, fct = modelList[[i]][["fct"]]),TRUE)
        }
        ),function(x) class(x)=="drc"))
      }
      
      non.convergence<-unique(unlist(drmModelListTmp))
      if(length(non.convergence)>0){
        bootData<-bootData[-non.convergence]
      }
      bootModelList <-list()
      for(i in 1:length(modelList)){
        bootModelList[[i]] <- lapply(bootData, function(x){
          suppressWarnings(
            eval(substitute(drm(formula, data = x, fct = fct0),
                            list(formula = modelList[[i]]$call$formula,
                                 fct0 = modelList[[i]][["fct"]]))
            ) # Fitting models using substitute is necessary for Stacking weights
          )
        }
        )
      }
      
      bootModelListTrans <- lapply(1:length(bootModelList[[1]]), function(i) lapply(bootModelList, "[[", i))
      
      # Compute weights on boot samples
      if(identical(modelWeights,"AIC")){
        AICList <-suppressWarnings(lapply(bootData, function(x) sapply(modelList, function(y) {AIC(my.fun(x,y))})))
        AICtmp <- do.call(rbind,AICList)
        modelWeights0 <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
      } else if(identical(modelWeights,"BIC")){
        BICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
        BICtmp <- do.call(rbind,BICList)
        modelWeights0 <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
      } else if(identical(modelWeights, "Stack")){  
        StackList <- lapply(bootModelListTrans, getStackingWeights)
        modelWeights0 <- do.call(cbind, StackList)
      } else {
        modelWeights0 <- do.call(cbind,rep(list(modelWeights),R))
      }
      
      # Estimate BMD in all models on boot samples
      bootbmdList<-list()
      for(i in 1:length(modelList)){
        bootbmdList[[i]] <- lapply(bootModelList[[i]], function(bootMod){
          try(bmd(bootMod,
                               bmr, backgType = backgType, backg = backg, def = def, interval = interval, 
                               display=FALSE, level=level)$Results[1,1], silent = TRUE)
        }
        )
      }
      bootbmdErrorList <- list()
      for(i in 1:length(modelList)){
        bootbmdErrorList[[i]] <- which(!sapply(bootbmdList[[i]],function(x) class(x)=="numeric"))
      }
      
      bmd.non.convergence<-unique(unlist(bootbmdErrorList))
      if(length(bmd.non.convergence) == length(bootData)){ 
        maBMDL <- NA 
        maBMDU <- NA
      } else {
        if(length(bmd.non.convergence)>0){
          for(i in 1:length(modelList)){
            bootbmdList[[i]]<-bootbmdList[[i]][-bmd.non.convergence]
            modelWeights0 <- modelWeights0[,-bmd.non.convergence]
          }
        }
        
        boot<-diag(t(matrix(unlist(bootbmdList), ncol = R - length(bmd.non.convergence), byrow = TRUE)) %*% modelWeights0)
        boot0<-boot[!is.na(boot)]
        
        if(bootInterval %in% c("percentile","Percentile")){
          maBMDL <- quantile(boot0,p=c(1-level), na.rm = TRUE) # ABC percentile lims.  
          maBMDU <- quantile(boot0,p=c(level), na.rm = TRUE)
        }
        if(identical(bootInterval,"BCa")){
          jackData <- list()
          for(i in 1:(dim(modelList[[1]]$data)[1])){
            jackData[[i]] <- modelList[[1]]$data[-i,]
          }
          bootJackList <-list()
          for(i in 1:length(modelList)){
            bootJackList[[i]] <- sapply(jackData, function(x){
              suppressWarnings(bmd(drm(modelList[[i]]$call$formula, data = x, type = modelList[[i]]$type, 
                                       fct = modelList[[i]][["fct"]]),
                                   bmr, backgType = backgType, backg = backg, def = def, interval = interval, 
                                   display=FALSE, level=level)$Results[1])
            }
            )
          }
          if(identical(modelWeights,"AIC")){
            AICJackList <-list()
            for(i in 1:length(modelList)){
              AICJackList[[i]] <- sapply(jackData, function(x){
                suppressWarnings(AIC(drm(modelList[[i]]$call$formula, data = x, type = modelList[[i]]$type, fct = modelList[[i]][["fct"]])))
              }
              )
            }
            AICtmp <- do.call(rbind,AICJackList)
            modelWeightsJack <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
          } else if(identical(modelWeights,"BIC")){
            BICJackList <-list()
            for(i in 1:length(modelList)){
              BICJackList[[i]] <- sapply(jackData, function(x){
                suppressWarnings(BIC(drm(modelList[[i]]$call$formula, data = x, type = modelList[[i]]$type, fct = modelList[[i]][["fct"]])))
              }
              )
            }
            BICtmp <- do.call(rbind,BICJackList)
            modelWeightsJack <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
          } else {
            modelWeightsJack <- do.call(cbind,rep(list(modelWeights),dim(modelList[[1]]$data)[1]))
          }
          
          bootjack<-diag(t(do.call(rbind,bootJackList)) %*% modelWeightsJack)
          
          maBMDL <- as.numeric(BCa(obs = maBMD, data = modelList[[1]]$data, boot0, bootjack)[1])
          maBMDU <- "Not available for BCa bootstrap"
        }
      }
    }
    if(identical(type,"curve")){
      bmrScaled0<-sum(sapply(bmdList, function(x){x$bmrScaled})*modelWeights0)
      maBMD <- bmdMACurve(modelList,modelWeights0,bmrScaled0)$Results[1]
      
      if(identical(bootstrapType,"nonparametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric")
      } else if(identical(bootstrapType,"semiparametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric")
      } else if(identical(bootstrapType,"parametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric")
      }
      
      drmModelListTmp <-list()
      for(i in 1:length(modelList)){
        drmModelListTmp[[i]] <- which(!sapply(lapply(bootData, function(x){
          try(drm(modelList[[i]]$call$formula, data = x, fct = modelList[[i]][["fct"]]),TRUE)
        }
        ),function(x) class(x)=="drc"))
      }
      
      non.convergence<-unique(unlist(drmModelListTmp))
      if(length(non.convergence)>0){
        bootData<-bootData[-non.convergence]
      }
      bootModelList <-list()
      for(i in 1:length(modelList)){
        bootModelList[[i]] <- lapply(bootData, function(x){
          suppressWarnings(
            eval(substitute(drm(formula, data = x, fct = fct0),
                            list(formula = modelList[[i]]$call$formula,
                                 fct0 = modelList[[i]][["fct"]]))
            ) # Fitting models using substitute is necessary for Stacking weights
          )
        }
        )
      }
      
      bootModelListTrans <- lapply(1:length(bootModelList[[1]]), function(i) lapply(bootModelList, "[[", i))
      
      if(identical(modelWeights,"AIC")){
        AICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(AIC(my.fun(x,y)))))
        AICtmp <- do.call(rbind,AICList)
        modelWeights0 <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
      } else if(identical(modelWeights,"BIC")){
        BICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
        BICtmp <- do.call(rbind,BICList)
        modelWeights0 <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
      } else if(identical(modelWeights, "Stack")){  
        modelWeightsList <- lapply(bootModelListTrans, getStackingWeights)
      } else {
        modelWeights0 <- do.call(cbind,rep(list(modelWeights),R))
      }
      
      if(!identical(modelWeights, "Stack")){ 
        modelWeightsList <- lapply(1:ncol(modelWeights0),function(i) modelWeights0[,i])
      }
      
      bootbmrList<-list()
      for(i in 1:length(modelList)){
        bootbmrList[[i]] <- sapply(bootModelList[[i]], function(bootMod){
          as.numeric(try(bmd(bootMod,
                               bmr, backgType = backgType, backg = backg, def = def, interval = interval, 
                               display=FALSE, level=level)$bmrScaled, silent = TRUE))
        }
        )
      }
      bootbmrListTrans <- lapply(1:length(bootbmrList[[1]]), function(i) sapply(bootbmrList, "[[", i))
      
      LLimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]]))[2]/10000
      ULimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]],decreasing=TRUE))[1]
      funk<-function(x,y,z){try(bmdMACurve(x,y,z,searchInterval=c(LLimit,ULimit))$Results[1],TRUE)}
      bmrScaledList<-as.list(rowSums(do.call(rbind,modelWeightsList)*do.call(rbind,bootbmrListTrans)))
      
      boot<-mapply(funk,bootModelListTrans,modelWeightsList,bmrScaledList)
      boot0<-suppressWarnings(as.numeric(boot[!is.na(as.numeric(boot))]))
      
      if(bootInterval %in% c("percentile","Percentile")){
        maBMDL <- quantile(boot0,p=c(1-level), na.rm = FALSE) # ABC percentile lims.  
        maBMDU <- quantile(boot0,p=c(level), na.rm = FALSE)
      }
      if(identical(bootInterval,"BCa")){
        jackData <- list()
        for(i in 1:(dim(modelList[[1]]$data)[1])){
          jackData[[i]] <- modelList[[1]]$data[-i,]
        }
        
        bootJackModelList <- lapply(jackData, function(x) lapply(modelList, function(y) suppressWarnings(my.fun(x,y))))
        
        
        if(identical(modelWeights,"AIC")){
          AICJackList <-lapply(jackData, function(x) sapply(modelList, function(y) suppressWarnings(AIC(my.fun(x,y)))))
          AICtmp <- do.call(rbind,AICJackList)
          modelWeightsJack <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
        } else if(identical(modelWeights,"BIC")){
          BICJackList <-lapply(jackData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
          BICtmp <- do.call(rbind,BICJackList)
          modelWeightsJack <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
        } else {
          modelWeightsJack <- do.call(cbind,rep(list(modelWeights),dim(modelList[[1]]$data)[1]))
        }
        modelWeightsJackList <- lapply(1:ncol(modelWeightsJack),function(i) modelWeightsJack[,i])
        
        jackbmrList<-list()
        for(i in 1:length(modelList)){
          jackbmrList[[i]] <- sapply(jackData, function(x){
            suppressWarnings(bmd(drm(modelList[[i]]$call$formula, data = x, type = modelList[[i]]$type, 
                                     fct = modelList[[i]][["fct"]]),
                                 bmr, backgType = backgType, backg = backg, def = def, interval = interval, 
                                 display=FALSE, level=level)$bmrScaled)
          }
          )
        }
        
        jackbmrListTrans <- lapply(1:length(jackbmrList[[1]]), function(i) sapply(jackbmrList, "[[", i))
        
        LLimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]]))[2]/10000
        ULimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]],decreasing=TRUE))[1]
        funk<-function(x,y,z){bmdMACurve(x,y,z,searchInterval=c(LLimit,ULimit))$Results[1]}
        bmrScaledJack<-as.list(rowSums(do.call(rbind,modelWeightsJackList)*do.call(rbind,jackbmrListTrans)))
        
        bootjack<-mapply(funk,bootJackModelList,modelWeightsJackList,bmrScaledJack)
        
        maBMDL <- as.numeric(BCa(obs = maBMD, data = modelList[[1]]$data, boot0, bootjack)[1])
        maBMDU <- "Not available for BCa bootstrap"
      }
    }
  }
  
  if(identical(modelList[[1]]$type,"binomial")){
    my.fun<-function(x,y){drm(y$call$formula, data = x, type="binomial", fct = y[["fct"]])}
    my.fun2<-function(x,y){drm(number~dose, data = x, type="binomial", fct = y[["fct"]])}
    
    
    if(identical(modelWeights,"AIC")){
      modelWeights0<-exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC))))/
        sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))))
    } else if(identical(modelWeights,"BIC")){
      modelWeights0<-exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC))))/
        sum(exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))))
    } else if(identical(modelWeights,"Stack")){
      modelWeights0<-getStackingWeights(modelList)
    } else {
      modelWeights0 <- modelWeights
    }
    
    if(identical(type,"Kang")){
      maBMD <- sum(modelWeights0 * sapply(bmdList, function(x){x$Results[,1]}))
      maBMDL <- sum(modelWeights0 * sapply(bmdList, function(x){x$Results[,2]}))
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
      if(identical(bootstrapType,"nonparametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric",aggregated = FALSE)
      } else if(identical(bootstrapType,"semiparametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric",aggregated = FALSE)
      } else if(identical(bootstrapType,"parametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric",aggregated = FALSE)
      }
      maBMD <- sum(modelWeights0 * sapply(bmdList, function(x){x$Results[,1]}))
      
      drmModelListTmp <-list()
      for(i in 1:length(modelList)){
        drmModelListTmp[[i]] <- which(!sapply(lapply(bootData, function(x){
          try(drm(modelList[[i]]$call$formula, data = x, type="binomial", fct = modelList[[i]][["fct"]]),TRUE)
        }
        ),function(x) class(x)=="drc"))
      }
      
      non.convergence<-unique(unlist(drmModelListTmp))
      if(length(non.convergence)>0){
        bootData<-bootData[-non.convergence]
      }
      bootModelList <-list()
      for(i in 1:length(modelList)){
        bootModelList[[i]] <- lapply(bootData, function(x){
          suppressWarnings(
            eval(substitute(drm(formula, data = x, fct = fct0, type = "binomial"),
                            list(formula = modelList[[i]]$call$formula,
                                 fct0 = modelList[[i]][["fct"]]))
            ) # Fitting models using substitute is necessary for Stacking weights
          )
        }
        )
      }
      
      bootModelListTrans <- lapply(1:length(bootModelList[[1]]), function(i) lapply(bootModelList, "[[", i))
      
      # Compute weights
      if(identical(modelWeights,"AIC")){
        AICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(AIC(my.fun(x,y)))))
        AICtmp <- do.call(rbind,AICList)
        modelWeights0 <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
      } else if(identical(modelWeights,"BIC")){
        BICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
        BICtmp <- do.call(rbind,BICList)
        modelWeights0 <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
      } else if(identical(modelWeights, "Stack")){  
        StackList <- lapply(bootModelListTrans, getStackingWeights)
        modelWeights0 <- do.call(cbind, StackList)
      } else {
        modelWeights0 <- do.call(cbind,rep(list(modelWeights),R))
      }
      
      # Estimate BMD in all models on boot samples
      bootbmdList<-list()
      for(i in 1:length(modelList)){
        bootbmdList[[i]] <- lapply(bootModelList[[i]], function(bootMod){
          try(bmd(bootMod,
                               bmr, backgType = backgType, backg = backg, def = def, interval = interval, 
                               display=FALSE, level=level)$Results[1,1], silent = TRUE)
        }
        )
      }
      bootbmdErrorList <- list()
      for(i in 1:length(modelList)){
        bootbmdErrorList[[i]] <- which(!sapply(bootbmdList[[i]],function(x) class(x)=="numeric"))
      }
      
      bmd.non.convergence<-unique(unlist(bootbmdErrorList))
      if(length(bmd.non.convergence) == length(bootData)){ 
        maBMDL <- NA 
        maBMDU <- NA
      } else {
        if(length(bmd.non.convergence)>0){
          for(i in 1:length(modelList)){
            bootbmdList[[i]]<-bootbmdList[[i]][-bmd.non.convergence]
          }
        }
        
        boot<-diag(t(matrix(unlist(bootbmdList), ncol = R - length(bmd.non.convergence), byrow = TRUE)) %*% modelWeights0[,-bmd.non.convergence])
        boot0<-boot[!is.na(boot)]
        
        if(bootInterval %in% c("percentile","Percentile")){
          maBMDL <- quantile(boot0,p=c(1-level), na.rm = FALSE) # ABC percentile lims.  
          maBMDU <- quantile(boot0,p=c(level), na.rm = FALSE)
        }
        if(identical(bootInterval,"BCa")){
          data.str <- modelList[[1]]$data
          data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
          data.e<-expandBinomial(data.str, 
                                 number = "number",
                                 total = "weights",
                                 dose = as.character(modelList[[1]]$call$formula[[3]]))
          jackData <- list()
          for(i in 1:(dim(data.e)[1])){
            jackData[[i]] <- data.e[-i,]
          }
          bootJackList <-list()
          for(i in 1:length(modelList)){
            bootJackList[[i]] <- sapply(jackData, function(x){
              suppressWarnings(bmd(my.fun2(x,modelList[[i]]),
                                   bmr, backgType = backgType, backg = backg, def = def, interval = interval, 
                                   display=FALSE, level=level)$Results[1])
            }
            )
          }
          if(identical(modelWeights,"AIC")){
            AICJackList <-list()
            for(i in 1:length(modelList)){
              AICJackList[[i]] <- sapply(jackData, function(x){
                suppressWarnings(AIC(my.fun2(x,modelList[[i]])))
              }
              )
            }
            AICtmp <- do.call(rbind,AICJackList)
            modelWeightsJack <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
          } else if(identical(modelWeights,"BIC")){
            BICJackList <-list()
            for(i in 1:length(modelList)){
              BICJackList[[i]] <- sapply(jackData, function(x){
                suppressWarnings(BIC(my.fun2(x,modelList[[i]])))
              }
              )
            }
            BICtmp <- do.call(rbind,BICJackList)
            modelWeightsJack <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
          } else {
            modelWeightsJack <- do.call(cbind,rep(list(modelWeights),dim(modelList[[1]]$data)[1]))
          }
          
          bootjack<-diag(t(do.call(rbind,bootJackList)) %*% modelWeightsJack)
          
          maBMDL <- as.numeric(BCa(obs = maBMD, data = data.e, boot0, bootjack)[1])
          maBMDU <- "Not available for BCa bootstrap"
        }
      }
    }
    if(identical(type,"curve")){
      bmrScaled0<-sum(sapply(bmdList, function(x){x$bmrScaled})*modelWeights0)
      maBMD <- bmdMACurve(modelList,modelWeights0,bmrScaled0)$Results[1]
      
      if(identical(bootstrapType,"nonparametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="nonparametric",aggregated = FALSE)
      } else if(identical(bootstrapType,"semiparametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="semiparametric",aggregated=FALSE)
      } else if(identical(bootstrapType,"parametric")){
        bootData <- bootDataGen(modelList[[1]],R=R,boot="parametric",aggregated=FALSE)
      }
      
      
      drmModelListTmp <-list()
      for(i in 1:length(modelList)){
        drmModelListTmp[[i]] <- which(!sapply(lapply(bootData, function(x){
          try(drm(modelList[[i]]$call$formula, data = x, type="binomial", fct = modelList[[i]][["fct"]]),TRUE)
        }
        ),function(x) class(x)=="drc"))
      }
      
      non.convergence<-unique(unlist(drmModelListTmp))
      if(length(non.convergence)>0){
        bootData<-bootData[-non.convergence]
      }
      
      bootModelList <-list()
      for(i in 1:length(modelList)){
        bootModelList[[i]] <- lapply(bootData, function(x){
          suppressWarnings(
            eval(substitute(drm(formula, data = x, fct = fct0, type = "binomial"),
                            list(formula = modelList[[i]]$call$formula,
                                 fct0 = modelList[[i]][["fct"]]))
            )
          )
        }
        )
      }
      
      bootModelListTrans <- lapply(1:length(bootModelList[[1]]), function(i) lapply(bootModelList, "[[", i))
      
      if(identical(modelWeights,"AIC")){
        AICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(AIC(my.fun(x,y)))))
        AICtmp <- do.call(rbind,AICList)
        modelWeights0 <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
        
      } else if(identical(modelWeights,"BIC")){
        BICList <-lapply(bootData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
        BICtmp <- do.call(rbind,BICList)
        modelWeights0 <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
      } else if(identical(modelWeights, "Stack")){  
        modelWeightsList <- lapply(bootModelListTrans, getStackingWeights)
      } else {
        modelWeights0 <- do.call(cbind,rep(list(modelWeights),R))
      }
      
      if(!identical(modelWeights, "Stack")){ 
        modelWeightsList <- lapply(1:ncol(modelWeights0),function(i) modelWeights0[,i])
      }
      
      bootbmrList<-list()
      for(i in 1:length(modelList)){
        bootbmrList[[i]] <- lapply(bootModelList[[i]], function(bootMod){
          try(bmd(bootMod,
                               bmr, backgType = backgType, backg = backg, def = def, interval = interval, 
                               display=FALSE, level=level)$bmrScaled, silent = TRUE)
        }
        )
      }
      bootbmrListTrans <- lapply(1:length(bootbmrList[[1]]), function(i) sapply(bootbmrList, "[[", i))
      
      bootbmrErrorList <- list()
      for(i in 1:length(modelList)){
        bootbmrErrorList[[i]] <- which(!sapply(bootbmrList[[i]],function(x) class(x)=="numeric"))
      }
      
      bmr.non.convergence<-unique(unlist(bootbmrErrorList))
      if(length(bmr.non.convergence) > 0){
        bootbmrListTrans <- bootbmrListTrans[-bmr.non.convergence]
        modelWeightsList <- modelWeightsList[-bmr.non.convergence]
        bootModelListTrans <- bootModelListTrans[-bmr.non.convergence]
      }
      
      LLimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]]))[2]/10000
      ULimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]],decreasing=TRUE))[1]
      funk<-function(x,y,z){
        try(bmdMACurve(x,y,z,searchInterval=c(LLimit,ULimit))$Results[1], silent = TRUE)}
      bmrScaledList<-as.list(rowSums(do.call(rbind,modelWeightsList)*do.call(rbind,bootbmrListTrans)))
      
      boot<-mapply(funk,bootModelListTrans,modelWeightsList,bmrScaledList)
      boot0<-suppressWarnings(as.numeric(boot[!is.na(as.numeric(boot))]))
      
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
        jackData <- list()
        for(i in 1:(dim(data.e)[1])){
          jackData[[i]] <- data.e[-i,]
        }
        
        bootJackModelList <- lapply(jackData, function(x) lapply(modelList, function(y) suppressWarnings(my.fun(x,y))))
        
        
        if(identical(modelWeights,"AIC")){
          AICJackList <-lapply(jackData, function(x) sapply(modelList, function(y) suppressWarnings(AIC(my.fun(x,y)))))
          AICtmp <- do.call(rbind,AICJackList)
          modelWeightsJack <- t(t(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp)))))/colSums(exp(-t(AICtmp - do.call(pmin, as.data.frame(AICtmp))))))
        } else if(identical(modelWeights,"BIC")){
          BICJackList <-lapply(jackData, function(x) sapply(modelList, function(y) suppressWarnings(BIC(my.fun(x,y)))))
          BICtmp <- do.call(rbind,BICJackList)
          modelWeightsJack <- t(t(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp)))))/colSums(exp(-t(BICtmp - do.call(pmin, as.data.frame(BICtmp))))))
        } else {
          modelWeightsJack <- do.call(cbind,rep(list(modelWeights),dim(modelList[[1]]$data)[1]))
        }
        modelWeightsJackList <- lapply(1:ncol(modelWeightsJack),function(i) modelWeightsJack[,i])
        
        jackbmrList<-list()
        for(i in 1:length(modelList)){
          jackbmrList[[i]] <- sapply(jackData, function(x){
            suppressWarnings(bmd(my.fun(x,modelList[[i]]),
                                 bmr, backgType = backgType, backg = backg, def = def, interval = interval, 
                                 display=FALSE, level=level)$bmrScaled)
          }
          )
        }
        
        jackbmrListTrans <- lapply(1:length(jackbmrList[[1]]), function(i) sapply(jackbmrList, "[[", i))
        
        LLimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]]))[2]/10000
        ULimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]],decreasing=TRUE))[1]
        funk<-function(x,y,z){bmdMACurve(x,y,z,searchInterval=c(LLimit,ULimit))$Results[1]}
        bmrScaledJack<-as.list(rowSums(do.call(rbind,modelWeightsJackList)*do.call(rbind,jackbmrListTrans)))
        
        bootjack<-mapply(funk,bootJackModelList,modelWeightsJackList,bmrScaledJack)
        
        maBMDL <- as.numeric(BCa(obs = maBMD, data = data.e, boot0, bootjack)[1])
        maBMDU <- "Not available for BCa bootstrap"
      }
    }
  }
  resMat<-matrix(c(maBMD,maBMDL),1,2)
  colnames(resMat) <- c("BMD_MA", "BMDL_MA")
  rownames(resMat) <- c("")
  
  used.Boot<-ifelse(identical(type,"bootstrap")|identical(type,"Bootstrap")|identical(type,"curve"),
                    length(boot0),NA)
  resBMD<-list(Results = resMat,
               Boot.samples.used = used.Boot,
               Interval = c(maBMDL, maBMDU))
  if(display){
    print(resMat)
  }
  
  class(resBMD) <- "bmd"
  invisible(resBMD)
}




