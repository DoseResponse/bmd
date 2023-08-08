computeWeightsFromSplit <- function(trainData, validateData, modelList){
  # Initialise weights to zero
  tmpWeights <- numeric(length(modelList))
  
  # Fit models to training data
  tmpModelList <- lapply(modelList, function(model){
    try(
      eval(substitute(drm(model$call$formula, weights = weights0, data = trainData, type = model$type, fct = model$fct),
       list(weights0 = model$call$weights))),
        silent = TRUE) # weights = eval(parse(text=as.character(model$call$weights)))
  })
  
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
  predMatrix <- sapply(tmpModelList, function(model) model$curve[[1]](validateData[, model$dataList$names$dName])) 
  
  # Convex optimization over the weights using CVXR
  alphaHat <- CVXR::Variable(length(tmpModelList))
  if(modelList[[1]]$type == "continuous"){
    objective <- CVXR::Minimize(sum((predMatrix %*% alphaHat - validateData[,modelList[[1]]$dataList$names$orName])^2))
  } else if (modelList[[1]]$type == "binomial"){
    objective <- CVXR::Minimize(sum((predMatrix %*% alphaHat - validateData[,as.character(modelList[[1]]$call$formula[[2]][[2]])])^2))
  }
  problem <- CVXR::Problem(objective, constraints = list(alphaHat <= 1, alphaHat >= 0,sum(alphaHat) == 1))
  result <- CVXR::solve(problem)
  res <- result$getValue(alphaHat)
  
  # Put optimised weights in correct place in weight vector
  tmpWeights[!convError] <- res
  tmpWeights <- tmpWeights * (tmpWeights > 0) # Slightly negative values can occur from CVXR
  tmpWeights <- tmpWeights / sum(tmpWeights) # Rescale weights
  tmpWeights
}

getDataSplits <- function(object){
  if(object$type == "continuous"){
    rowNum <- 1:object$sumList$lenData
    split1_indices <- aggregate(
      rowNum, list(object$dataList$dose), 
      FUN = function(rowNum){
        # if uneven number of repetitions, choose randomly if floor(rep/2) or ceiling(rep/2) is chosen.
        tmpSize <- sample(c(floor(length(rowNum)/2),ceiling(length(rowNum)/2)), 1) 
        sample(rowNum, size = tmpSize)
      })$x |> unlist()
    
    splitList <- list(split1 = object$origData[split1_indices,],
                      split2 = object$origData[-split1_indices,])
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

getSuperLearnerWeights <- function(modelList){
  get_w_error <- TRUE
  tries <- 1
  while(get_w_error && (tries<10)){
    dataSplits <- getDataSplits(modelList[[1]])
    
    weightsT1V2 <- computeWeightsFromSplit(dataSplits[[1]], dataSplits[[2]], modelList) |> try(silent = TRUE)
    weightsT2V1 <- computeWeightsFromSplit(dataSplits[[2]], dataSplits[[1]], modelList) |> try(silent = TRUE)
    get_w_error <- any(c(inherits(weightsT1V2, "try-error"), inherits(weightsT2V1, "try-error")))
    tries <- tries + 1
  }
  if(get_w_error){return(rep(NA, length(modelList)))}
  else{
    returnWeights <- drop(weightsT1V2 + weightsT2V1) / 2
    returnWeights
  }
}
