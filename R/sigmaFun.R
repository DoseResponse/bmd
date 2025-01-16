sigmaFun <- function(object, formula){
  # Assertions
  if(!class(object) == "drc"){
    stop('object must be a dose-response model of class "drc" ')
  }
  if(length(unique(object$dataList$curveid)) != 1){
    stop("dose-response models with multiple curves not supported for heteroscedasticity analysis")
  }
  
  if(class(formula) != "formula"){
    stop('argument "formula" must be of class "formula"')
  }
  
  # Add fitted values and residuals to data
  data <- object$data |>
    dplyr::mutate(fitted = fitted(object),
                  residuals = residuals(object))
  
  # Aggregate data
  data.agg <- data |>
    dplyr::group_by(fitted) |>
    dplyr::summarise(dose0 = mean(.data[[object$dataList$names$dName]]),
                     sigma0 = sqrt(mean(residuals^2)))
  colnames(data.agg)[2] <- object$dataList$names$dName
  
  
  formula <- as.formula(formula)
  formula0 <- reformulate(attr(terms(formula), "term.labels"), response = "sigma0")
  
  sigma.mod <- lm(formula0, data = data.agg)
  
  ret.fun <- function(x){
    newdata0 <- data.frame(dose0 = x, fitted = object$curve[[1]](x))
    colnames(newdata0)[1] <- object$dataList$names$dName
    
    predict(sigma.mod, newdata0)
  }
  
  # Checking for roots.
  # NOT STABLE IF THERE ARE MULTIPLE ROOTS IN DOSE RANGE!
  interval0 <- range(data[[object$dataList$names$dName]], na.rm = TRUE)
  root.try <- try(uniroot(ret.fun, 
                          interval = interval0), silent = TRUE)
  
  if(!inherits(root.try, "try-error")){
    stop("Root detected in variance function. Choose a different model for the variance. \n")
  }
  
  # Return object
  ret.list <- list(ret.fun = ret.fun, sigma.mod = sigma.mod, data.agg = data.agg, model = object)
  class(ret.list) <- "drc.sigma.fun"
  ret.list
}

plot.drc.sigma.fun <- function(object, gridsize = 300){
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
  polygonY <- c(curveFun(dosePts) + 1.96*object$ret.fun(dosePts), 
                rev(curveFun(dosePts) - 1.96*object$ret.fun(dosePts)) )
  
  p1 <- ggplot() +
    geom_polygon(aes(x = polygonX, y = polygonY), alpha = 0.1) +
    geom_line(aes(x = dosePts, y = curveFun(dosePts))) +
    geom_point(aes(x = dose, y = resp)) +
    scale_x_continuous(trans = "pseudo_log") +
    labs(x = doseName, y = respName)
  
  p2 <- ggplot(object$data.agg) +
    geom_point(aes(x = .data[[dName]], y = sigma0)) +
    geom_function(fun = object$ret.fun) +
    scale_x_continuous(trans = "pseudo_log")
  
  (gridExtra::grid.arrange(p1, p2))
  invisible(list(p1,p2))
}