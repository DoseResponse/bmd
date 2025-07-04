#' Dose response modeling with heterogeneous variance
#' 
#' Fit a dose-response model with heterogeneous variance dependending on dose
#' level.
#' 
#' The aim to provide an R package calculating the benchmark dose (BMD) and the
#' lower limit of the corresponding 95\% confidence interval (BMDL) for
#' continuous and quantal dose-response data for a range of dose-response
#' models based on the available definitions of the benchmark dose concepts.
#' 
#' REFERENCES TO BE ADDED/WRITTEN
#' 
#' @param formula formula for the dose-response relationship
#' @param var.formula one-sided formula specifying the dependance of the dose
#' values and/or the fitted values on the point-wise standard error
#' @param data data.frame containing the observations
#' @param fct a list with three or more elements specifying the non-linear
#' function, the accompanying self starter function, the names of the parameter
#' in the non-linear function and, optionally, the first and second derivatives
#' as well as information used for calculation of ED values. Currently
#' available functions include, among others, the four- and five-parameter
#' log-logistic models LL.4, LL.5 and the Weibull model W1.4. Use
#' drc::getMeanFunctions for a full list.
#' @param curveStart numerical of length equal to the number of parameters for
#' the curve. Starting values for the curve parameters (optional).
#' @return dose-response model with a heterogeneous variance structure of class
#' \code{drcHetVar}.
#' 
#' The primary objective is to use this model for benchmark dose estimation
#' based on the hybrid method with a heterogeneous variance structure. This can
#' be done using the \code{bmdHetVar} function.
#' 
#' A plot method is available, which can be useful for assessing the fit of the
#' variance structure.
#' @author Signe M. Jensen and Jens Riis Baalkilde
#' @keywords models nonlinear
#' @examples
#' 
#' library(drc)
#' library(drcData)
#' library(bmd)
#' # install.packages("gridExtra") # OPTIONAL - USED FOR PLOTTING A drcHetVar OBJECT.
#' 
#' # ryegrass data
#' set.seed(123)
#' ryegrass.LL.4.hetVar <- drmHetVar(rootl ~ conc, ~ fitted + I(fitted^2),
#'                                   data = ryegrass, fct = LL.4())
#' plot(ryegrass.LL.4.hetVar)
#' bmdHetVar(ryegrass.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1,
#'           def = "hybridExc", R = 50, level = 0.95, progressInfo = TRUE, display = TRUE)
#' bmdHetVar(ryegrass.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, 
#'           def = "hybridExc", R = 50, level = 0.95, 
#'           bootType = "parametric", progressInfo = TRUE, display = TRUE) # parametric bootstrap
#' 
#' # barley data
#' set.seed(123)
#' barley.LL.4.hetVar <- drmHetVar(weight ~ Dose, ~ fitted + I(fitted^2), data = barley, fct = LL.4())
#' plot(barley.LL.4.hetVar)
#' 
#' # GiantKelp data
#' set.seed(123)
#' GiantKelp.LL.4.hetVarSq <- drmHetVar(tubeLength ~ dose, ~ fitted + I(fitted^2), 
#'                                      data = GiantKelp, fct = LL.4())
#' plot(GiantKelp.LL.4.hetVarSq)
#' 
#' GiantKelp.LL.4.hetVarLogSq <- drmHetVar(tubeLength ~ dose, ~ log(dose+1) + I(log(dose+1)^2), 
#'                                         data = GiantKelp, fct = LL.4())
#' plot(GiantKelp.LL.4.hetVarLogSq)
#' 
#' 
#' 
#' @export
drmHetVar <- function(formula, var.formula, data, fct, curveStart = NULL) {
  call <- match.call()
  
  if(!inherits(formula, "formula")){
    stop('argument "formula" must be of class "formula"')
  }
  
  if(!inherits(var.formula, "formula")){
    stop('argument "var.formula" must be of class "formula"')
  }
  
  if(missing(data)){
    stop('argument "data" must be supplied')
  }
  
  if(missing(fct)){
    stop('argument "fct" must be supplied')
  }
  
  if(!requireNamespace("dplyr")){
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
  start <- unlist(drmHetVarSelfStarter(formula, var.formula, mf, fct, curveStart))
  
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
