#' Comparison of benchmark dose estimates
#' 
#' Comparison of benchmark doses and benchmark dose estimates in models with 
#' multiple curve levels.
#' 
#' DEVELOPMENT VERSION - USE AT OWN RISK AND PLEASE PROVIDE FEEDBACK ABOUT USAGE
#' 
#' 
#' @param object object of class \code{bmd}
#' @param operator additive differences ("-") or proportional differences ("/")
#' @param vcov. which function to use to compute variance-covariance matrix. Only
#' relevant when object is bmd() object, and not bmdBoot()
#' @param display whether to display the results
#' 
#' @return A matrix with columns containing the estimates, estimated standard 
#' errors, values of t-statistics and p-values for the null hypothesis that the 
#' ratio equals 1 or that the difference equals 0 (depending on the operator 
#' argument).
#' 
#' @author Signe M. Jensen and Jens Riis Baalkilde
#' 
#' @keywords models nonlinear
#' @examples
#' 
#' library(drc)
#' 
#' ## Examples to be added
#' 
#' @export
compParmBmd <- function(object, operator = c("-", "/"), vcov. = vcov, 
                         display = TRUE){
  # ASSERTIONS TO BE ADDED
  
  operator <- match.arg(operator)
  
  if(is.null(object$Boot.samples.used)){
    type <- "delta"
  } else {
    type <- "boot"
  }
  
  
  presentVec <- rownames(object$Results)
  lenPV <- length(presentVec)
  if (lenPV < 2) {
    stop("No parameters to compare")
  }
  
  lenRV <- lenPV * (lenPV - 1)/2
  cpMat <- matrix(0, lenRV, 4)
  compParm <- rep("", lenRV)
  
  if(identical(type, "delta")){
    bmr <- object$bmdOpts$bmr
    backgType <- object$bmdOpts$backgType
    backg <- object$bmdOpts$backg
    controlSD <- object$bmdOpts$controlSD
    def <- object$bmdOpts$def
    respTrans <- object$bmdOpts$respTrans
    
    EDlist <- bmd.edfct(object$model)
    parmMat <- object$model$parmMat
    nCurves <- lenPV # ncol(parmMat)
    
    # bmrScaledList
    bmrScaledList <- getBmrScaledList(object$model, bmr, backgType, backg, controlSD, def, respTrans)
    
    fctDerivx <- object$model$fct$derivx
    
    varMat <- vcov.(object$model)
    pmodelsMatrixList <- getpmodelsMatrixList(object$model)
    
    dBmdFun <- function(iCurve){
      parmChosen <- parmMat[, iCurve]
      
      EDeval <- EDlist(parmChosen, bmrScaledList$bmrScaledMat[iCurve,], type = "absolute", reference = "control") 
      bmdVal <- EDeval[[1]]
      dBmdVal <- (EDeval[[2]] + bmrScaledList$dBmrScaled[,iCurve] / fctDerivx(bmdVal, t(parmMat))[iCurve]) %*% pmodelsMatrixList[[iCurve]]
      dBmdVal
    }
    
    if (identical(operator, "/")) {
      hypVal <- 1
      fct <- function(val1, val2) {
        ifelse(val2 == 0, 
               NA,
               val1/val2) 
      }
      dfct <- function(ind1, ind2) {
        if(object$Results[presentVec[ind2],1] == 0){
          NA
        } else {
          dBmd1 <- dBmdFun(ind1)
          dBmd2 <- dBmdFun(ind2)
          
          # transVec <- c(1/parm[ind[2]], -parm[ind[1]]/(parm[ind[2]]^2))
          transVec <- dBmd1/object$Results[presentVec[ind1],1] -dBmd2 * object$Results[presentVec[ind1],1]/(object$Results[presentVec[ind2],1]^2)
          SEval <- sqrt(transVec %*% varMat %*% t(transVec))
          SEval
        }
      }
    }
    if (identical(operator, "-")) {
      hypVal <- 0
      fct <- function(val1, val2) {
        val1 - val2
      }
      # transVec <- c(1, -1)
      dfct <- function(ind1, ind2) {
        dBmd1 <- dBmdFun(ind1)
        dBmd2 <- dBmdFun(ind2)
        transVec <- dBmd1 - dBmd2
        SEVal <- sqrt(transVec %*% varMat %*% t(transVec))
        SEVal
      }
    }
    
    
    degfree <- df.residual(object$model)
    if (is.null(degfree)) {
      degfree <- 100
    }
    if (object$model$type == "continuous") {
      pFct <- function(x) {
        pt(x, degfree)
      }
    } else {
      pFct <- pnorm
    }
    
    k <- 1
    for (i in 1:lenPV) {
      for (j in 1:lenPV) {
        if (j <= i) {
          next
        }
        cpMat[k, 1] <- fct(object$Results[presentVec[i],1], object$Results[presentVec[j],1]) # fct(presentVec[c(i, j)])
        cpMat[k, 2] <- dfct(i,j) # dfct(presentVec[c(i, j)])
        tVal <- (cpMat[k, 1] - hypVal)/cpMat[k, 2]
        cpMat[k, 3] <- tVal
        cpMat[k, 4] <- pFct(-abs(tVal)) + (1 - pFct(abs(tVal)))
        compParm[k] <- paste("BMD:", paste(presentVec[c(i, j)], collapse = operator))
        k <- k + 1
      }
    }
  } else if(identical(type, "boot")){
    if (identical(operator, "/")) {
      hypVal <- 1
      fct <- function(val1, val2) {
        ifelse(val2 == 0, 
               NA,
               val1/val2) 
      }
    }
    if (identical(operator, "-")) {
      hypVal <- 0
      fct <- function(val1, val2) {
        val1 - val2
      }
    }
    
    k <- 1
    for (i in 1:lenPV) {
      for (j in 1:lenPV) {
        if (j <= i) {
          next
        }
        bootEst <- fct(object$bootEst[,presentVec[i]], object$bootEst[,presentVec[j]])
        cpMat[k, 1] <- fct(object$Results[presentVec[i],1], object$Results[presentVec[j],1])
        cpMat[k, 2] <- sd(bootEst, na.rm = TRUE)
        cpMat[k, 3] <- NA 
        lowTail <- mean(bootEst <= hypVal, na.rm = TRUE)
        uppTail <- mean(bootEst >= hypVal, na.rm = TRUE)
        cpMat[k, 4] <- 2 * min(lowTail, uppTail)
        compParm[k] <- paste("BMD:", paste(presentVec[c(i, j)], collapse = operator))
        k <- k + 1
      }
    }
  }
  
  
  dimnames(cpMat) <- list(compParm, c("Estimate", "Std. Error", 
                                      "t-value", "p-value"))
  if (display) {
    cat("\nComparison of BMD", "\n\n")
    printCoefmat(cpMat)
  }
  invisible(cpMat)
}