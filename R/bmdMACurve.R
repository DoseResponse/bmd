#' Calculate Model-Averaged BMD Curve, helper to bmdMA
#'
#' @description
#' Helper function for `bmdMA` that computes Benchmark Dose (BMD) values using 
#' model averaging. The function creates a weighted average of multiple dose-response 
#' models and finds the dose corresponding to a specified benchmark response level.
#'
#' @param modelList list; A list of fitted dose-response models from which to compute 
#'        the model average. Each model should be a fitted object with components 
#'        including `fct`, `parmMat`, `dataList`, etc.
#' @param modelWeights numeric vector; Weights for model averaging, typically derived 
#'        from model selection criteria (e.g., AIC weights). Must sum to 1 and have 
#'        the same length as `modelList`
#' @param bmrScaled0 numeric; The benchmark response level(s) for which to calculate 
#'        BMD. For single curves, a scalar value; for multiple curves, a vector with 
#'        length equal to the number of curves
#' @param searchInterval character or numeric; Search interval for root finding. 
#'        If "dataBased" (default), uses data range with small lower bound adjustment. 
#'        If numeric, should be a vector of length 2 for single curves, or a matrix 
#'        with 2 columns for multiple curves specifying `[lower, upper]` bounds
#'
#' @details
#' The function handles two scenarios:
#' 
#' **Single Curve Analysis:**
#' Creates a weighted combination of model functions:
#' \deqn{f_{MA}(x) = \sum_{i=1}^{m} w_i \cdot f_i(x)}{f_MA(x) = sum of w_i * f_i(x)}
#' where \eqn{w_i}{w_i} are model weights and \eqn{f_i(x)}{f_i(x)} are individual model functions.
#' 
#' **Multiple Curve Analysis:**
#' Computes BMD separately for each curve using the same model averaging approach 
#' but applied to curve-specific data and parameters.
#' 
#' The BMD is found by solving:
#' \deqn{f_{MA}(BMD) = BMR}{f_MA(BMD) = BMR}
#' using root finding algorithms.
#' 
#' **Model Function Construction:**
#' The function dynamically constructs model functions by:
#' \itemize{
#'   \item Extracting function templates from `mtList()`
#'   \item Substituting fixed parameters and estimated coefficients
#'   \item Handling special cases like Fractional Polynomial models
#'   \item Creating weighted combinations of all models
#' }
#'
#' @return An object of class "bmd" containing:
#' \itemize{
#'   \item \strong{Results}: A matrix with BMD values. For single curves, a 1x1 matrix; 
#'         for multiple curves, an nx1 matrix where n is the number of curves
#'   \item \strong{MACurve}: The model-averaged function used for BMD calculation. 
#'         For single curves, returns the combined function; for multiple curves, returns NULL
#' }
#'
#' @section Search Interval:
#' When `searchInterval = "dataBased"`:
#' \itemize{
#'   \item Lower bound: Second smallest dose value divided by 10,000
#'   \item Upper bound: Maximum dose value in the dataset
#' }
#' 
#' Custom intervals can be provided to avoid convergence issues or focus on 
#' specific dose ranges of interest.
#'
#' @section Model Support:
#' The function supports various dose-response models including:
#' \itemize{
#'   \item Standard 4-parameter and 5-parameter models
#'   \item Fractional Polynomial models (special handling)
#'   \item Models with fixed parameters
#' }
#'
#' @note
#' This is primarily an internal helper function for `bmdMA`. Direct usage requires 
#' properly formatted model lists and weights. The function assumes all models in 
#' `modelList` are compatible and fitted to the same dataset structure.
#'
#' @seealso 
#' \code{\link{bmdMA}} for the main model averaging function,
#' \code{\link{uniroot}} for the root finding algorithm used internally
#'
#' @examples
#' \dontrun{
#' # Typically called internally by bmdMA, but can be used directly:
#' 
#' # Assume you have a list of fitted models and weights
#' models <- list(model1, model2, model3)  # fitted drc models
#' weights <- c(0.5, 0.3, 0.2)  # AIC weights
#' bmr <- 0.1  # 10% benchmark response
#' 
#' # Calculate model-averaged BMD
#' result <- bmdMACurve(models, weights, bmr)
#' 
#' # Extract BMD value
#' bmd_value <- result$Results[1, 1]
#' 
#' # Use the averaged curve function
#' curve_function <- result$MACurve
#' }
#'
#' @export
bmdMACurve<-function(modelList,modelWeights,bmrScaled0, searchInterval="dataBased"){
  nCurves <- ncol(modelList[[1]]$parmMat)
  
  if(nCurves == 1){
    fList <- mtList()
    
    f.all<-list()
    fi<-lapply(modelList,function(x){fList[[x[["fct"]][["name"]]]]$fct})
    
    for(i in 1:length(modelList)){
      if(!identical(modelList[[i]]$fct$text,"Fractional polynomial")){
      if(length(modelList[[i]]$fct$fixed)==5){
        parm<-modelList[[i]]$fct$fixed
        parm[is.na(parm)]<-coef(modelList[[i]])
        body(fi[[i]])<-gsub("b1",parm[1],paste(body(fi[[i]])))[2]
        body(fi[[i]])<-gsub("c1",parm[2],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("d1",parm[3],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("e1",parm[4],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("f1",parm[5],paste(body(fi[[i]])))
        f.all[[i]] <- fi[[i]]
      } else if(length(modelList[[i]]$fct$fixed)==4){
        parm<-modelList[[i]]$fct$fixed
        parm[is.na(parm)]<-coef(modelList[[i]])
        body(fi[[i]])<-gsub("b1",parm[1],paste(body(fi[[i]])))[2]
        body(fi[[i]])<-gsub("c1",parm[2],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("d1",parm[3],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("e1",parm[4],paste(body(fi[[i]])))
        f.all[[i]] <- fi[[i]]
      } 
    } else if(identical(modelList[[i]]$fct$text,"Fractional polynomial")){
        fi[[i]] <- fList[["FPL.4"]]$fct
        parm<-modelList[[i]]$fct$fixed
        parm[is.na(parm)]<-coef(modelList[[i]])
        body(fi[[i]])<-gsub("b1",parm[1],paste(body(fi[[i]])))[2]
        body(fi[[i]])<-gsub("c1",parm[2],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("d1",parm[3],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("e1",parm[4],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("p1",unlist(strsplit(modelList[[i]]$fct$name, "\\,|\\(|\\)"))[2],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("p2",unlist(strsplit(modelList[[i]]$fct$name, "\\,|\\(|\\)"))[3],paste(body(fi[[i]])))
        f.all[[i]] <- fi[[i]]
      }
    }
    
    Body<-paste(mapply(FUN=function(x,y){paste(y, " * (",body(x),")")}, f.all, modelWeights), collapse = " + ")
    
    args <- as.character("dose")
    
    eval(parse(text = paste('g <- function(', args, ') { return(' , Body , "-", bmrScaled0, ')}', sep='')))
    
    if(identical(searchInterval,"dataBased")){
    LLimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]]))[2]/10000
    ULimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]],decreasing=TRUE))[1]
    } else { 
      LLimit <- searchInterval[1]
      ULimit <- searchInterval[2]
    }
    BMD<-uniroot(g,interval=c(LLimit,ULimit))$root
    
    resMat<-matrix(c(BMD),1,1)
    colnames(resMat) <- c("BMD")
    rownames(resMat) <- c("")
  } else {
    curveNames <- colnames(modelList[[1]]$parmMat)
    uniqueCurveids <- unique(modelList[[1]]$dataList$curveid)
    nMods <- length(modelList)
    
    curveFcts <- lapply(modelList, 
                        function(model){ 
                          function(x) model$curve[[1]](x)
                        })
    
    getBMDi <- function(i){
      colIndex <- which(curveNames[i] == uniqueCurveids)
      
      g <- function(x){
        curveVal <- colSums(modelWeights * matrix(sapply(1:nMods, function(j) curveFcts[[j]](x)[, colIndex]), nrow = nMods, byrow = TRUE))
        curveVal - bmrScaled0[i]
      }
      
      if(identical(searchInterval,"dataBased")){
        LLimit <- unique(sort(modelList[[1]]$dataList$dose[modelList[[1]]$dataList$curveid == curveNames[i]]))[2]/10000
        ULimit <- unique(sort(modelList[[1]]$dataList$dose[modelList[[1]]$dataList$curveid == curveNames[i]],decreasing=TRUE))[1]
      } else { 
        LLimit <- searchInterval[i,1]
        ULimit <- searchInterval[i,2]
      }
      
      BMD <- uniroot(g,interval=c(LLimit,ULimit))$root
      BMD
    }
    
    BMD <- sapply(1:nCurves, getBMDi)
    
    resMat <- matrix(BMD, nrow = nCurves, ncol = 1)
    colnames(resMat) <- "BMD"
    rownames(resMat) <- colnames(modelList[[1]]$parmMat)
    g <- NULL
  }
  
  resBMD<-list(Results = resMat,
               MACurve = g)
  class(resBMD) <- "bmd"
  invisible(resBMD)
}


