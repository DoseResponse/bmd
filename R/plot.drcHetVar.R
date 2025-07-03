
#' @title S3 method
#' @export
plot.drcHetVar <- function(x, ...){
  object <- x
  dots <- list(...)
  if (!is.null(dots$gridsize)){
    gridsize <- dots$gridsize
  } else {
    gridsize <- 300
  }
  
  # Add assertion of gridExtra
  if(!requireNamespace("gridExtra")){
    stop('package "gridExtra" must be installed to plot drcHetVar object')
  }
  
  if(!requireNamespace("ggplot2")){
    stop('package "ggplot2" must be installed to plot drcHetVar object')
  }
  
  if(!requireNamespace("dplyr")){
    stop('package "dplyr" must be installed to plot drcHetVar object')
  }
  
  # Plot of model
  dose <- object$dataList[["dose"]]
  resp <- object$dataList[["resp"]]
  doseName <- object$dataList$names$dName
  respName <- object$dataList$names$rName
  
  xLimits <- range(dose)
  xLimits0 <- pmax(xLimits, 1e-8)
  dosePts <- c(0,exp(seq(log(xLimits0[1]), log(xLimits0[2]), length = gridsize-1)))
  dosePts[1] <- max(xLimits[1],0)
  dosePts[gridsize] <- xLimits[2]    
  
  curveFun <- object$curve
  
  polygonX <- c(dosePts, rev(dosePts))
  polygonY <- c(curveFun(dosePts) + 1.96*object$sigmaFun(dosePts), 
                rev(curveFun(dosePts) - 1.96*object$sigmaFun(dosePts)) )
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_polygon(aes(x = polygonX, y = polygonY), alpha = 0.1) +
    ggplot2::geom_line(aes(x = dosePts, y = curveFun(dosePts))) +
    ggplot2::geom_point(aes(x = dose, y = resp)) +
    ggplot2::scale_x_continuous(trans = "pseudo_log") +
    ggplot2::labs(x = doseName, y = respName)
  
  df <- data.frame(dose = dose, resp = resp, residuals = object$residuals) # fitted.values = object$fitted.values, 
  df.agg <- dplyr::summarise(dplyr::group_by(df, dose), 
                               sigma0 = sqrt(mean(residuals^2)), groups = ".keep")
  
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_point(aes(x = dose, y = sigma0), data = df.agg) +
    ggplot2::geom_line(aes(x = dosePts, y = object$sigmaFun(dosePts))) +
    ggplot2::scale_x_continuous(trans = "pseudo_log") +
    ggplot2::labs(x = doseName, y = "sqrt[mean(residuals^2)]")
  
  (gridExtra::grid.arrange(p1, p2))
  invisible(list(p1,p2))
}
