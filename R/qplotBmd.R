qplotBmd <- function(x, ..., interval = c("BMDL", "twosided", "none"), col = FALSE, add = FALSE){
  object <- x
  
  if(!inherits(object, "bmd")){
    stop('qplotBmd only works for plotting objects of type "bmd"')
  }
  
  if(!is.null(object$modelWeights)){
    stop('qplotBmd does not for for model-averaged BMD')
  }
  
  
  model <- object$model
  interval <- match.arg(interval)
  
  if(nrow(object$Results) == 1){
    # One curve
    xVert <- rep(object$Results[1],2)
    yVert <- c(model$curve[[1]](object$Results[1]), 0)
    
    xHoriz <- c(0, ifelse(interval != "none", object$interval[1], object$Results[1]))
    yHoriz <- rep(model$curve[[1]](object$Results[1]),2)
    
    xInt <- c(object$interval[1], ifelse(interval == "twosided", object$interval[2], object$Results[1]))
    xLow <- rep(object$interval[1],2)
    xUpp <- rep(object$interval[2],2)
    
    returnLayers <- list(geom_line(aes(x = xHoriz, y = yHoriz), linetype = 2),
                         geom_line(aes(x = xVert, y = yVert)))
    
    if(interval == "BMDL"){
      returnLayers <- c(returnLayers,
                        geom_line(aes(x = xInt, y = yHoriz)),
                        geom_line(aes(x = xLow, y = yVert)))
    } else if(interval == "twosided"){
      returnLayers <- c(returnLayers,
                        geom_line(aes(x = xInt, y = yHoriz)),
                        geom_line(aes(x = xLow, y = yVert)),
                        geom_line(aes(x = xUpp, y = yVert)))
    }
  } else {
    # Multiple curves
    curveLevels <- unique(model$dataList$curveid)
    
    xVert <- rep(object$Results[,1],each = 2)
    curveID <- names(xVert)
    yVertVal <- diag(model$curve[[1]](object$Results[,1])[,match(rownames(object$Results), curveLevels)])
    yVert <- numeric(nrow(object$Results)*2)
    yVert[1:nrow(object$Results)*2 - 1] <- yVertVal
    
    xHoriz <- unlist(lapply(1:nrow(object$Results), function(row) c(0, ifelse(interval != "none", object$interval[row, 1], object$Results[row, 1]))))
    yHoriz <- rep(diag(model$curve[[1]](object$Results[,1])[,match(rownames(object$Results), curveLevels)]), each = 2)
    
    xInt <- unlist(lapply(1:nrow(object$Results), function(row) c(object$interval[row, 1], ifelse(interval == "twosided", object$interval[row, 2], object$Results[row, 1]))))
    xLow <- rep(object$interval[,1],each = 2)
    xUpp <- rep(object$interval[,2], each = 2)
    
    if(col){
      returnLayers <- list(geom_line(aes(x = xHoriz, y = yHoriz, col = curveID, group = curveID), linetype = 2),
                           geom_line(aes(x = xVert, y = yVert, col = curveID, group = curveID)))
    } else {
      returnLayers <- list(geom_line(aes(x = xHoriz, y = yHoriz, linetype = curveID, group = curveID)),
                           geom_line(aes(x = xVert, y = yVert, linetype = curveID, group = curveID)))
    }
    
    if(interval == "BMDL" & col){
      returnLayers <- c(returnLayers,
                        geom_line(aes(x = xInt, y = yHoriz, col = curveID, group = curveID)),
                        geom_line(aes(x = xLow, y = yVert, col = curveID, group = curveID)))
    } else if(interval == "BMDL" & !col){
      returnLayers <- c(returnLayers,
                        geom_line(aes(x = xInt, y = yHoriz, linetype = curveID, group = curveID)),
                        geom_line(aes(x = xLow, y = yVert, linetype = curveID, group = curveID)))
    } else if(interval == "twosided" & col){
      returnLayers <- c(returnLayers,
                        geom_line(aes(x = xInt, y = yHoriz, col = curveID, group = curveID)),
                        geom_line(aes(x = xLow, y = yVert, col = curveID, group = curveID)),
                        geom_line(aes(x = xUpp, y = yVert, col = curveID, group = curveID)))
    } else if(interval == "twosided" & !col){
      returnLayers <- c(returnLayers,
                        geom_line(aes(x = xInt, y = yHoriz, linetype = curveID, group = curveID)),
                        geom_line(aes(x = xLow, y = yVert, linetype = curveID, group = curveID)),
                        geom_line(aes(x = xUpp, y = yVert, linetype = curveID, group = curveID)))
    }
  }
  if(add){
    returnLayers
  } else {
    qplotDrc(model, ..., col = col) +
      returnLayers
  }
}
