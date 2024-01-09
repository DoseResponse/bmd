MACurve <- function(x, modelList, modelWeights){
  # compute weights
  if(identical(modelWeights,"AIC")){
    modelWeights0<-exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC))))/
      sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))))
  } else if(identical(modelWeights,"BIC")){
    modelWeights0<-exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC))))/
      sum(exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))))
  } else {
    modelWeights0 <- modelWeights
  }
  
  sapply(x, function(x0){
    vals <- sapply(modelList, function(mod) mod$curve[[1]](x0))
    sum(modelWeights0 * vals)
  })
}
