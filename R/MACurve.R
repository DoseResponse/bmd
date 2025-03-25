MACurve <- function(x, modelList, modelWeights, stackingSeed = 1, stackingSplits = 2){
  
  
  # compute weights
  if(identical(modelWeights,"AIC")){
    modelWeights0<-exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))/2)/
      sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))/2))
  } else if(identical(modelWeights,"BIC")){
    modelWeights0<-exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))/2)/
      sum(exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))/2))
  } else if(identical(modelWeights, "Stack")){
    # If stackingSeed supplied, save initial seed for later, and set seed for stacking
    if (!is.null(stackingSeed)) {
      sysSeed <- .GlobalEnv$.Random.seed
      set.seed(stackingSeed, kind = "Mersenne-Twister", normal.kind = "Inversion")
    }
    
    # estimate weights
    modelWeights0 <- getStackingWeights(modelList, nSplits = stackingSplits)
    
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
  
  sapply(x, function(x0){
    vals <- sapply(modelList, function(mod) mod$curve[[1]](x0))
    sum(modelWeights0 * vals)
  })
}
