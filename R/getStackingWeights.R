computeWeightsFromSplit <- function(trainData, validateData, modelList){
  
  # Fit models to training data
  if(is.null(modelList[[1]]$call$pmodels)){
    tmpModelList <- lapply(modelList, function(model){
      try(#update(model, data = trainData),
        eval(substitute(drm(model$call$formula, weights = weights0, curveid = curveid0,
                            data = trainData, type = model$type, fct = model$fct, control = drmc(noMessage = TRUE)),
         list(weights0 = model$call$weights,
              curveid0 = model$call$curveid)
         )),
          silent = TRUE)
    })
  } else {
    tmpModelList <- lapply(modelList, function(model){
      try(#update(model, data = trainData),
        eval(substitute(drm(model$call$formula, weights = weights0, curveid = curveid0,pmodels = pmodels0,
                            data = trainData, type = model$type, fct = model$fct, control = drmc(noMessage = TRUE)),
                        list(weights0 = model$call$weights,
                             curveid0 = model$call$curveid,
                             pmodels0 = model$call$pmodels)
        )),
        silent = TRUE)
    })
  }
  
  # Only compute weights for succesfully converged models
  convError <- sapply(tmpModelList, function(mod_try) inherits(mod_try, "try-error"))
  tmpModelList <- tmpModelList[!convError]
  
  # Treat data differently for binomial data
  if(modelList[[1]]$type == "binomial"){
    data.e<-expandBinomial(validateData, 
                           number = as.character(modelList[[1]]$call$formula[[2]][[2]]),
                           total = as.character(modelList[[1]]$call$formula[[2]][[3]]),
                           dose = as.character(modelList[[1]]$call$formula[[3]]))
    validateData <- data.e
  }
  
  if(ncol(modelList[[1]]$parmMat) == 1){
    predMatrix <- sapply(tmpModelList, function(model) model$curve[[1]](validateData[[model$dataList$names$dName]]))
    predMatrix <- matrix(predMatrix, ncol = length(tmpModelList))
  } else {
    predVec <- function(model, dose, curveid){
      mat <- model$curve[[1]](dose)
      sapply(1:length(dose), 
             function(i){
               val <- mat[i, which(curveid[i] == unique(model$dataList$curveid))]
               if(is.na(val)){ val <- mean(mat[i,], na.rm = TRUE)}
               val
               })
    }
    predMatrix <- sapply(tmpModelList, 
                       function(model){ predVec(model, 
                                                validateData[[model$dataList$names$dName]], 
                                                validateData[[model$call$curveid]])
                         }
    )
  }
  
  # Convex optimization over the weights using CVXR
  alphaHat <- CVXR::Variable(length(tmpModelList))
  if(modelList[[1]]$type == "continuous"){
    objective <- CVXR::Minimize(sum((predMatrix %*% alphaHat - matrix(validateData[[modelList[[1]]$dataList$names$orName]]))^2))
  } else if (modelList[[1]]$type == "binomial"){
    objective <- CVXR::Minimize(sum((predMatrix %*% alphaHat - validateData[[as.character(modelList[[1]]$call$formula[[2]][[2]])]])^2))
  }
  problem <- CVXR::Problem(objective, constraints = list(alphaHat <= 1, alphaHat >= 0,sum(alphaHat) == 1))
  result <- CVXR::solve(problem)
  res <- result$getValue(alphaHat)
  
  # Initialise weights to zero
  tmpWeights <- numeric(length(modelList))
  # Put optimised weights in correct place in weight vector
  tmpWeights[!convError] <- res
  tmpWeights <- tmpWeights * (tmpWeights > 0) # Slightly negative values can occur from CVXR
  tmpWeights <- tmpWeights / sum(tmpWeights) # Rescale weights
  tmpWeights
}

getDataSplits <- function(object, nSplits){
  if(object$type == "continuous"){
    rowNum <- 1:object$sumList$lenData
    # split1_indices <- aggregate(
    #   rowNum, list(object$dataList$dose, object$dataList$curveid), 
    #   FUN = function(rowNum){
    #     # if uneven number of repetitions, choose randomly if floor(rep/2) or ceiling(rep/2) is chosen.
    #     tmpSize <- sample(c(floor(length(rowNum)/2),ceiling(length(rowNum)/2)), 1) 
    #     sample(rowNum, size = tmpSize)
    #   })$x |> unlist()
    # 
    # splitList <- list(split1 = object$origData[split1_indices,],
    #                   split2 = object$origData[-split1_indices,])
    
    # splitNum <- aggregate(rowNum, list(object$dataList$dose, object$dataList$curveid), 
    #                       function(x) rep(sample(1:nSplits), ceiling(length(x)/nSplits))[1:length(x)]
    #                       )$x |> unlist()
    
    splitAvailable <- table(rep(sample(1:nSplits), ceiling(object$sumList$lenData/nSplits))[1:object$sumList$lenData])
    
    curveLevels <- sample(unique(object$dataList$curveid))
    doseLevels <- sample(unique(object$dataList$dose))
    
    splitNum <- integer(object$sumList$lenData)
    
    for(iCurve in 1:length(curveLevels)){
      for(jDose in 1:length(doseLevels)){
        rowNumDoseCurve <- rowNum[(object$dataList$dose == doseLevels[jDose]) & (object$dataList$curveid == curveLevels[iCurve])]
        splitsAssigned <- integer(0)
        missingVals <- length(rowNumDoseCurve)
        while(missingVals > 0){
          splitsToSample <- names(splitAvailable[splitAvailable == max(unique(splitAvailable))])
          addVals <- sample(splitsToSample, size = min(missingVals, length(splitsToSample)))
          splitAvailable[names(table(addVals))] <- splitAvailable[names(table(addVals))] - table(addVals)
          splitsAssigned <- c(splitsAssigned, as.numeric(addVals))
          missingVals <- length(rowNumDoseCurve) - length(splitsAssigned)
        }
        
        splitNum[rowNumDoseCurve] <- splitsAssigned
      }
    }
    splitList <- split(object$origData, splitNum)
  } else if(object$type == "binomial"){
    data.str <- object$data
    data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
    data.e<-expandBinomial(data.str, 
                           number = "number",
                           total = "weights",
                           dose = as.character(object$call$formula[[3]]))
    data.e$rowNum <-1:nrow(data.e)
    
    split1_indices <- aggregate(
      data.e$rowNum, list(data.e[,as.character(object$call$formula[[3]])]),
      FUN = function(rowNum){
        # if uneven number of repetitions, choose randomly if floor(rep/2) or ceiling(rep/2) is chosen.
        tmpSize <- sample(c(floor(length(rowNum)/2),ceiling(length(rowNum)/2)), 1) 
        sample(rowNum, size = tmpSize)
      })$x |> unlist()
    
    data.e.splitList <- list(data.e.split1 = data.e[split1_indices,],
                      data.e.split2 = data.e[-split1_indices,])
    
    splitList <- lapply(data.e.splitList,
                            function(data){
                              tmp <- aggregate(data[,c("number", "weights")], 
                                               by = list((data[,as.character(object$call$formula[[3]])])), 
                                               sum)
                              colnames(tmp) <- c(as.character(object$call$formula[[3]]),
                                                 as.character(object$call$formula[[2]])[[2]],
                                                 as.character(object$call$formula[[2]])[[3]])
                              tmp
                            })
  }
  
  splitList
}

getStackingWeights <- function(modelList, nSplits = 2){
  get_w_error <- TRUE
  tries <- 1
  while(get_w_error && (tries<10)){
    dataSplits <- getDataSplits(modelList[[1]], nSplits)
    
    # weightsT1V2 <- computeWeightsFromSplit(dataSplits[[1]], dataSplits[[2]], modelList) |> try(silent = TRUE)
    # weightsT2V1 <- computeWeightsFromSplit(dataSplits[[2]], dataSplits[[1]], modelList) |> try(silent = TRUE)
    # get_w_error <- any(c(inherits(weightsT1V2, "try-error"), inherits(weightsT2V1, "try-error")))
    weightList <- lapply(1:nSplits, 
                         function(i) {
                           computeWeightsFromSplit(do.call(rbind, dataSplits[-i]), 
                                                   dataSplits[[i]], modelList) |> try(silent = TRUE)
                           })
    
    get_w_error <- any(sapply(weightList, function(x) inherits(x, "try-error")))
    tries <- tries + 1
  }
  if(get_w_error){return(rep(NA, length(modelList)))}
  else{
    returnWeights <- colSums(do.call(rbind, weightList)) / nSplits  # drop(weightsT1V2 + weightsT2V1) / 2
    returnWeights
  }
}
