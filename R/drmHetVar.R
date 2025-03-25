drmHetVar <- function(object, var.formula){
  # Assertions
  if(!class(object) == "drc"){
    stop('object must be a dose-response model of class "drc" ')
  }
  if(length(unique(object$dataList$curveid)) != 1){
    stop("dose-response models with multiple curves not supported for heteroscedasticity analysis")
  }
  
  if(class(var.formula) != "formula"){
    stop('argument "formula" must be of class "formula"')
  }
  
  if(!require("dplyr")){
    stop('package "dplyr" must be installed to fit dose-response model with heterogeneous variance')
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
  
  
  formula <- as.formula(var.formula)
  formula0 <- reformulate(attr(terms(formula), "term.labels"), response = "sigma0")
  
  sigma.mod <- lm(formula0, data = data.agg)
  
  sigma.fun <- function(x){
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
  ret.list <- list(model = object, sigmaFun = sigma.fun, var.formula = var.formula, sigmaMod = sigma.mod, data.agg = data.agg)
  class(ret.list) <- "drcHetVar"
  ret.list
}

