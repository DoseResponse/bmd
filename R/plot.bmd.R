#' Plotting benchmark doses
#' 
#' \code{plot.bmd} displays benchmark dose values with options to plot
#' confidence intervals as well.
#' 
#' This function is a simple function to plot benchmark dose values along with
#' the fitted curve.
#' 
#' @aliases plot.bmd \method{plot}{bmd}
#' @param x an object of class 'bmd'.
#' @param ... arguments to be passed on to \code{plot.drc}, if \code{add =
#' FALSE}
#' @param interval option to plot only the lower limit of the confidence
#' interval for the benchmark dose ("BMDL", default), both limit of the
#' confidence interval ("twosided"), or no confidence interval ("none").
#' @return Creates a plot. No value returned.
#' @author Jens Riis Baalkilde.
#' @keywords ggplot
#' @examples
#' 
#' library(drc)
#' library(drcData)
#' 
#' ## Fitting model and calculating BMD. 
#' model <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#' bmd0 <- bmd(model, bmr = 0.1, backgType = "modelBased", def = "relative")
#' 
#' ## Plotting
#' plot(bmd0)
#' 
#' # Plot both limits of confidence interval
#' plot(bmd0, interval = "twosided")
#' 
#' # Pass argument to plot.bmd to plot confidence band around curve
#' plot(bmd0, type = "confidence")
#' 
#' 
#' 
plot.bmd <- function(x, ..., interval = c("BMDL", "twosided", "none")){
  object <- x
  model <- object$model
  interval <- match.arg(interval)
  
  p0 <- plot(model, ...)
  
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
