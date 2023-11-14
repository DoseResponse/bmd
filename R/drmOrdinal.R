drmOrdinal <- function(levels, dose, weights, data, fct, p.epsilon = 10^(-16)){
  # list of merged levels
  levelsMerged <- list()
  for(i in 2:(length(levels))){
    levelsMerged[[i-1]] <- paste(levels[i:length(levels)], collapse = '+')
  } 
  
  # fit model for each element of levelsMerged
  drmList <- list()
  drmList <- lapply(levelsMerged, function(cat){ 
    formula.string <- paste0("(",cat,")/",weights, "~", dose)
    eval(substitute(drm(formula0, weights = eval(parse(text=weights)), data = data, type = "binomial", fct = fct),
                    list(formula0 = eval(parse(text=formula.string)))))
  }
  )
  names(drmList) <- unlist(levelsMerged)
  
  pFun <- function(x){
    prob <- sapply(levels, function(cat){
      cat.i <- which(cat == levels)
      if(cat.i == 1){
        pmax(p.epsilon, 1 - drmList[[1]]$curve[[1]](x))
      } else if(cat.i == length(levels)){
        pmax(p.epsilon, drmList[[cat.i-1]]$curve[[1]](x))
      } else{
        pmax(p.epsilon, drmList[[cat.i-1]]$curve[[1]](x) - drmList[[cat.i]]$curve[[1]](x))
      }
    }
    )
    names(prob) <- levels
    prob
  }
  
    
  res.list <- list(drmList, levels, levelsMerged, dose, weights, data, fct, pFun)
  names(res.list) <- c("drmList", "levels", "levelsMerged", "dose", "weights", "data", "fct", "pFun")
  class(res.list) <- "drcOrdinal"
  
  res.list
}
