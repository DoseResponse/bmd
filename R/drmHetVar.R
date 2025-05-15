drmHetVar <- function(formula, var.formula, data, fct) {
  call <- match.call()
  
  if(class(formula) != "formula"){
    stop('argument "formula" must be of class "formula"')
  }
  
  if(class(var.formula) != "formula"){
    stop('argument "var.formula" must be of class "formula"')
  }
  
  if(missing(data)){
    stop('argument "data" must be supplied')
  }
  
  if(missing(fct)){
    stop('argument "fct" must be supplied')
  }
  
  if(!require("dplyr")){
    stop('package "dplyr" must be installed to fit dose-response model with heterogeneous variance')
  }
  
  # Extract response and dose
  mf <- model.frame(formula, data)
  mf <- mf[complete.cases(mf),]
  rName <- as.character(formula)[2]
  dName <- as.character(formula)[3]
  resp <- mf[[rName]]
  dose <- mf[[dName]]
  
  dataList <- list(dose = dose, resp = resp, names = list(dName = dName, rName = rName))
  
  # Define curve function (mean)
  curveFun <- function(x, par) {
    fct$fct(x, t(par))
  }
  
  # Define tau function constructor
  makeTauFun <- function(var.formula) {
    function(x, tauPar, curvePar) {
      fitted <- curveFun(x, curvePar)
      env <- data.frame(x = x, "fitted" = fitted)
      colnames(env)[1] <- dName
      design_matrix <- model.matrix(var.formula, data = env)
      as.vector(design_matrix %*% tauPar)
    }
  }
  
  # Build tau function
  tauFun <- makeTauFun(var.formula)
  
  # Total number of parameters:
  # - Mean model (from fct, usually fixed length)
  # - Variance model (length depends on model matrix)
  # tmp_design_matrix <- model.matrix(var.formula, data = data.frame(dName = dose, "fitted" = resp))
  # n_var_par <- ncol(tmp_design_matrix)
  n_mean_par <- sum(is.na(fct$fixed)) # length(fct$names)
  
  # Negative log-likelihood
  negLogLik <- function(par) {
    curvePar <- par[1:n_mean_par]
    tauPar <- par[-(1:n_mean_par)] 
    
    mu <- curveFun(dose, curvePar)
    sigma <- tauFun(dose, tauPar, curvePar)
    sigmaSq <- sigma^2
    
    if (any(sigma <= 0)) return(1e10)  # enforce positivity
    
    sum(log(sigmaSq)) + sum( (resp - mu)^2 / sigmaSq)
  }
  
  # Initial values (somewhat naive)
  # start_curve <- coef(drm(formula, data = data, fct = fct)) # rep(mean(y), n_mean_par)
  # start_tau <- rep(sd(y), n_var_par)
  # start <- c(start_curve, 0.05)# start_tau)
  start <- unlist(drmHetVarSelfStarter(formula, var.formula, data, fct))
  
  fit <- optim(start, negLogLik, method = "BFGS", hessian = TRUE)
  
  # Unpack results
  curvePar <- fit$par[1:n_mean_par]
  names(curvePar) <- fct$names
  sigmaPar <- fit$par[-(1:n_mean_par)]
  tmp_env <- data.frame(dName = dose, "fitted" = resp)
  colnames(tmp_env)[1] <- dName
  names(sigmaPar) <- colnames(model.matrix(var.formula, data = tmp_env))
  
  # Define sigmaFun
  curve <- function(x) curveFun(x, par = curvePar)
  sigmaFun <- function(x) tauFun(x, tauPar = sigmaPar, curvePar = curvePar)
  
  # sumList
  sumList <- list(numObs = length(resp),
                  numPar = length(curvePar) + length(sigmaPar))
  
  # Output
  object <- list(
    curvePar = curvePar,
    sigmaPar = sigmaPar,
    value = fit$value,
    convergence = fit$convergence,
    message = fit$message,
    hessian = fit$hessian,
    curve = curve,
    sigmaFun = sigmaFun,
    formula = formula,
    var.formula = var.formula,
    fct = fct,
    data = data,
    dataList = dataList,
    sumList = sumList,
    call = call,
    fitted.values = curveFun(dose, fit$par[1:n_mean_par]),
    residuals = resp - curveFun(dose, fit$par[1:n_mean_par])
  )
  
  class(object) <- c("drcHetVar", "drc")
  
  return(object)
}