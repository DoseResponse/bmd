getpmodelsMatrixList <- function(object){
  curveLevels <- colnames(object$parmMat)
  numLevels <- ncol(object$parmMat)
  parNames <- object$parNames
  numCurvePar <- length(unique(parNames[[2]]))
  numPar <- length(coef(object))
  curveParNames <- parNames[[2]]
  uniqCurveParNames <- unique(parNames[[2]])
  parSuffix <- parNames[[3]]
  
  
  pmodelsMatrixList <- list()
  if (numLevels == 1) 
  {
    mat0 <- diag(rep(1, numPar))
    rownames(mat0) <- uniqCurveParNames
    colnames(mat0) <- parNames[[1]]
    pmodelsMatrixList[[1]] <- mat0
  } else {
    indexFun <- function(row,col){
      (uniqCurveParNames[row] == curveParNames[col]) * 
        ((parSuffix[col] == "(Intercept)") + grepl(curveLevels[i], parSuffix[col]))
    }
    for(i in 1:numLevels){
      mat0 <- outer(1:numCurvePar, 1:numPar, indexFun)
      rownames(mat0) <- uniqCurveParNames
      colnames(mat0) <- parNames[[1]]
      
      pmodelsMatrixList[[i]] <- mat0
    }
    names(pmodelsMatrixList) <- curveLevels
  }
  
  pmodelsMatrixList
}
