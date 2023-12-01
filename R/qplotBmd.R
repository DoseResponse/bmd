qplotBmd <- function(x, ..., interval = c("BMDL", "twosided", "none"), add = FALSE){
  object <- x
  interval <- match.arg(interval)
  
  xVert <- rep(object$Results[1],2)
  yVert <- c(model$curve[[1]](object$Results[1]), 0)
  
  xHoriz <- c(0, ifelse(interval != "none", object$interval[1], object$Results[1]))
  yHoriz <- rep(model$curve[[1]](object$Results[1]),2)
  
  xInt <- c(object$interval[1], ifelse(interval == "twosided", object$interval[2], object$Results[1]))
  xLow <- rep(object$interval[1],2)
  xUpp <- rep(object$interval[2],2)
  
  returnLayers <- list(ggplot2:::geom_line(aes(x = xHoriz, y = yHoriz), linetype = 2),
                       ggplot2:::geom_line(aes(x = xVert, y = yVert)))
  
  if(interval == "BMDL"){
    returnLayers <- c(returnLayers,
                      ggplot2:::geom_line(aes(x = xInt, y = yHoriz)),
                      ggplot2:::geom_line(aes(x = xLow, y = yVert)))
  } else if(interval == "twosided"){
    returnLayers <- c(returnLayers,
                      ggplot2:::geom_line(aes(x = xInt, y = yHoriz)),
                      ggplot2:::geom_line(aes(x = xLow, y = yVert)),
                      ggplot2:::geom_line(aes(x = xUpp, y = yVert)))
  }
  
  if(add){
    returnLayers
  } else {
    qplotDrc(object$model, ...) +
      returnLayers
  }
}
