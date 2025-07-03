#' Fitting ordinal dose-response models
#' 
#' A model fitting function for analysis of ordinal dose-response data.
#' 
#' This functions fits a dose-response model for ordinal dose-response data by
#' fitting a series of binomial dose-response models.
#' 
#' @param levels a character vector of the levels of the response variable in
#' increasing order, as they appear in the supplied data set
#' @param dose a character string specifying the column with the dose values in
#' the supplied data set
#' @param weights a character string specifying the column containing the
#' number of observations pr. group
#' @param blocks a character string specifying the column containing the blocks
#' of the experiment, if available (optional)
#' @param data a dataframe containing the observations
#' @param fct a list with three or more elements specifying the non-linear
#' function, the accompanying self starter function, the names of the parameter
#' in the non-linear function and, optionally, the first and second derivatives
#' as well, for the individual fitted curves. For more information see the help
#' page for the "drm" function in the "drc" package.  Currently available
#' functions for ordinal dose-response models include, among others, the
#' log-logistic models \code{\link[drc]{LL.4}}, the log-normal model
#' \code{\link[drc]{LN.4}} and the two Weibull models \code{\link[drc]{W1.4}}
#' and \code{\link[drc]{W2.4}}. Use \code{\link[drc]{getMeanFunctions}} for a
#' full list.
#' @param p.epsilon numeric value specifying the lower bound for the
#' probabilites for each level returned by the function pFun created when
#' fitting the model. Default value is 10^(-16)
#' @return An object of (S3) class 'drcOrdinal'.
#' @author Signe M. Jensen and Jens Riis Baalkilde
#' @references # ADD REFERENCES
#' @keywords models nonlinear
#' @examples
#' 
#' library(drc)
#' library(drcData)
#' data(guthion)
#' guthionS <- subset(guthion, trt == "S")
#' 
#' guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total",
#'                           dose = "dose", data = guthionS, fct = LL.2())
#' 
#' plot(guthionS.LL, xlim = c(15,55)) # uses ggplot
#' 
#' @export
drmOrdinal <- function(levels, dose, weights, blocks, data, fct, p.epsilon = 1e-16){
  if(missing(blocks)){
    blocks = NULL
  }
  # list of merged levels
  levelsMerged <- list()
  for(i in 2:(length(levels))){
    levelsMerged[[i-1]] <- paste(levels[i:length(levels)], collapse = '+')
  } 
  
  # fit model for each element of levelsMerged
  drmList <- list()
  drmList <- lapply(levelsMerged, function(cat){ 
    formula.string <- paste0("(",cat,")/",weights, "~", dose)
    eval(substitute(drm(formula0, weights = eval(parse(text=weights)), data = data, type = "binomial", fct = fct),
                    list(formula0 = eval(parse(text=formula.string)))))
  }
  )
  names(drmList) <- unlist(levelsMerged)
  
  pFun <- function(x){
    prob <- sapply(levels, function(cat){
      cat.i <- which(cat == levels)
      if(cat.i == 1){
        pmax(p.epsilon, 1 - drmList[[1]]$curve[[1]](x))
      } else if(cat.i == length(levels)){
        pmax(p.epsilon, drmList[[cat.i-1]]$curve[[1]](x))
      } else{
        pmax(p.epsilon, drmList[[cat.i-1]]$curve[[1]](x) - drmList[[cat.i]]$curve[[1]](x))
      }
    }
    )
    names(prob) <- levels
    prob
  }
  
    
  res.list <- list(drmList, levels, levelsMerged, dose, weights, blocks, data, fct, pFun)
  names(res.list) <- c("drmList", "levels", "levelsMerged", "dose", "weights", "blocks", "data", "fct", "pFun")
  class(res.list) <- "drcOrdinal"
  
  res.list
}
