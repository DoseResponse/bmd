computeWeightsFromSplit <- function(trainData, validateData, modelList){
  # Initialise weights to zero
  tmpWeights <- numeric(length(modelList))
  
  # Fit models to training data
  tmpModelList <- lapply(modelList, function(model){
    try(drm(model$call$formula, data = trainData, type = model$type, fct = model[["fct"]]),
        silent = TRUE)
  })
  
  # Only compute weights for succesfully converged models
  convError <- sapply(tmpModelList, function(mod_try) inherits(mod_try, "try-error"))
  tmpModelList <- tmpModelList[!convError]
  predMatrix <- sapply(tmpModelList, function(model) model$curve[[1]](validateData[, model$dataList$names$dName])) 
  
  # Convex optimization over the weights using CVXR
  alphaHat <- CVXR::Variable(length(tmpModelList))
  objective <- CVXR::Minimize(sum((predMatrix %*% alphaHat - validateData[,modelList[[1]]$dataList$names$orName])^2))
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
  rowNum <- 1:object$sumList$lenData
  split1_indices <- aggregate(
    rowNum, list(object$dataList$dose), 
    FUN = function(rowNum){
      # if uneven number of repetitions, choose randomly if floor(rep/2) or ceiling(rep/2) is chosen.
      tmpSize <- sample(c(floor(length(rowNum)/2),ceiling(length(rowNum)/2)), 1) 
      sample(rowNum, size = tmpSize)
    })$x |> unlist()
  
  list(split1 = object$origData[split1_indices,],
       split2 = object$origData[-split1_indices,])
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
