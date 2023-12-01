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
  
  if(add){
    returnLayers
  } else {
    ggplotDrc(object$model, ...) +
      returnLayers
  }
}
