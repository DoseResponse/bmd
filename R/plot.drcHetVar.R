plot.drcHetVar <- function(object, gridsize = 300){
  # Add assertion of gridExtra
  
  dName <- colnames(object$data.agg)[2]
  
  # Plot of model
  dose <- object$model$dataList[["dose"]]
  resp <- object$model$dataList[["origResp"]]
  doseName <- object$model$dataList$names$dName
  respName <- object$model$dataList$names$orName
  
  xLimits <- range(dose)
  xLimits0 <- pmax(xLimits, 1e-8)
  dosePts <- c(0,exp(seq(log(xLimits0[1]), log(xLimits0[2]), length = gridsize-1)))
  dosePts[1] <- max(xLimits[1],0)
  dosePts[gridsize] <- xLimits[2]    
  
  curveFun <- object$model$curve[[1]]
  
  polygonX <- c(dosePts, rev(dosePts))
  polygonY <- c(curveFun(dosePts) + 1.96*object$sigmaFun(dosePts), 
                rev(curveFun(dosePts) - 1.96*object$sigmaFun(dosePts)) )
  
  p1 <- ggplot() +
    geom_polygon(aes(x = polygonX, y = polygonY), alpha = 0.1) +
    geom_line(aes(x = dosePts, y = curveFun(dosePts))) +
    geom_point(aes(x = dose, y = resp)) +
    scale_x_continuous(trans = "pseudo_log") +
    labs(x = doseName, y = respName)
  
  p2 <- ggplot(object$data.agg) +
    geom_point(aes(x = .data[[dName]], y = sigma0)) +
    geom_function(fun = object$sigmaFun) +
    scale_x_continuous(trans = "pseudo_log")
  
  (gridExtra::grid.arrange(p1, p2))
  invisible(list(p1,p2))
}