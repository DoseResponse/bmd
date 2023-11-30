plot.bmd <- function(x, ..., interval = c("BMDL", "twosided", "none")){
  object <- x
  interval <- match.arg(interval)
  
  p0 <- plot(object$model, ...)
  
  xVert <- rep(object$Results[1],2)
  yVert <- c(model$curve[[1]](object$Results[1]), 0)
  
  xHoriz <- c(min(p0[,1]), ifelse(interval != "none", object$interval[1], object$Results[1]))
  yHoriz <- rep(model$curve[[1]](object$Results[1]),2)
  
  xInt <- c(object$interval[1], ifelse(interval == "twosided", object$interval[2], object$Results[1]))
  xLow <- rep(object$interval[1],2)
  xUpp <- rep(object$interval[2],2)
  
  lines(xHoriz, yHoriz, lty = 2)
  lines(xVert, yVert)
  
  if(interval == "BMDL"){
    lines(xInt, yHoriz)
    lines(xLow, yVert)
  } else if(interval == "twosided"){
    lines(xInt, yHoriz)
    lines(xLow, yVert)
    lines(xUpp, yVert)
  }
}
