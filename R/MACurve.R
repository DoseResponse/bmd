#' Model-average dose-response curves
#' 
#' Computing weighted average response estimates across multiple dose-response
#' curves.
#' 
#' The aim to provide an R package calculating the benchmark dose (BMD) and the
#' lower limit of the corresponding 95\% confidence interval (BMDL) for
#' continuous and quantal dose-response data for a range of dose-response
#' models based on the available definitions of the benchmark dose concepts.
#' 
#' Details on the implemented definitions and methods can be found in Crump
#' (2002)
#' 
#' @param x a vector of dose values for which the weighted average of response
#' estimates are to be computed
#' @param modelList list of models of class \code{drc}
#' @param modelWeights character string specifying the type of weights used,
#' "AIC", "BIC" or "Stack", or a vector of the same length as the modelList
#' with user defined weights
#' @param stackingSeed integer or NULL: Random seed to use in the data split in
#' the estimation of the Stacking Weights, when \code{modelWeights = "Stack"}.
#' The global seed is reset to the initial value after estimation of the
#' weights, so this option does not interfere with a globally set seed.
#' @param stackingSplits integer or "LOO": When \code{modelWeights = "Stack"},
#' the Stacking Weights are estimated, which are based on V-fold
#' cross-validation. The stackingSplits argument sets the number V of data
#' splits used in the cross validation. The "LOO" (Leave one out) is a shortcut
#' to setting V equal to the number of observations in the data set.
#' @return numeric
#' @author Jens Riis Baalkilde
#' @keywords models nonlinear model averaging
#' @examples
#' 
#' library(bmd)
#' library(drc)
#' library(drcData)
#' library(ggplot2)
#' 
#' # fit models to aconiazide data
#' aconiazide.LL.3 <- drm(weightChange ~ dose,data = aconiazide,fct = LL.3())
#' aconiazide.LN.3 <- drm(weightChange ~ dose,data = aconiazide,fct = LN.3())
#' aconiazide.W1.3 <- drm(weightChange ~ dose,data= aconiazide,fct = W1.3())
#' aconiazide.W2.3 <- drm(weightChange ~ dose,data= aconiazide,fct = W2.3())
#' 
#' # plot the MA curve
#' plot(aconiazide.LL.3, type = "obs")
#' curve(
#'   MACurve(x, modelList = list(aconiazide.LL.3, aconiazide.LN.3,aconiazide.W1.3, aconiazide.W2.3),
#'           modelWeights = "AIC"),
#'   add = TRUE)
#' 
#' # or plot using ggplot2
#' qplotDrc(aconiazide.LL.3, type = "obs") +
#'  geom_function(fun = function(x){ 
#'                 MACurve(x, modelList = list(aconiazide.LL.3, aconiazide.LN.3,
#'                                             aconiazide.W1.3, aconiazide.W2.3), 
#'                         modelWeights = "AIC")
#'                 })
#' 
#' @export
MACurve <- function(x, modelList, modelWeights, stackingSeed = 1, stackingSplits = 2){
  
  
  # compute weights
  if(identical(modelWeights,"AIC")){
    modelWeights0<-exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))/2)/
      sum(exp(-(sapply(modelList,AIC)-min(sapply(modelList,AIC)))/2))
  } else if(identical(modelWeights,"BIC")){
    modelWeights0<-exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))/2)/
      sum(exp(-(sapply(modelList,BIC)-min(sapply(modelList,BIC)))/2))
  } else if(identical(modelWeights, "Stack")){
    # If stackingSeed supplied, save initial seed for later, and set seed for stacking
    if (!is.null(stackingSeed)) {
      sysSeed <- .GlobalEnv$.Random.seed
      set.seed(stackingSeed, kind = "Mersenne-Twister", normal.kind = "Inversion")
    }
    
    # estimate weights
    modelWeights0 <- getStackingWeights(modelList, nSplits = stackingSplits)
    
    # If stackingSeed supplied, restore initial seed
    if (!is.null(stackingSeed)) {
      if (!is.null(sysSeed)) {
        .GlobalEnv$.Random.seed <- sysSeed 
      } else {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }
  } else {
    modelWeights0 <- modelWeights
  }
  
  sapply(x, function(x0){
    vals <- sapply(modelList, function(mod) mod$curve[[1]](x0))
    sum(modelWeights0 * vals)
  })
}
