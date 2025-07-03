#' Plotting benchmark doses using ggplot2
#' 
#' \code{qplotBmd} displays benchmark dose values with options to plot
#' confidence intervals as well using \code{ggplot2}.
#' 
#' This function is a simple function to plot benchmark dose values
#' 
#' @param x an object of class 'bmd'.
#' @param ... arguments to be passed on to qplotDrc, if \code{add = FALSE}
#' @param interval option to plot only the lower limit of the confidence
#' interval for the benchmark dose ("BMDL", default), both limit of the
#' confidence interval ("twosided"), or no confidence interval ("none").
#' @param col logical. If TRUE then multiple curves specified by "curveid" in
#' the dose-response model are distinguised by colours rather than point shapes
#' and line types
#' @param add logical. If TRUE then the functions returns a list of plot layers
#' to be added to an already existing ggplot.
#' @return A \code{ggplot} object. If the option \code{add} is used, a list of
#' \code{ggplot} layers is returned.
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
#' # Plot
#' qplotBmd(bmd0, interval = "twosided", add = FALSE)
#' 
#' qplotDrc(model,type="confidence") +
#'   qplotBmd(bmd0, interval = "twosided", add = TRUE)
#' 
#' qplotBmd(bmd0, interval = "twosided", add = FALSE)
#' 
#' 
#' 
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
