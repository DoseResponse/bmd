bmdProfileCIgrid <- function(object, bmr, backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"), backg, controlSD, def = c("excess", "additional", "relative", "extra", "added", "hybridExc", "hybridAdd", "point"),
                         level = 0.95, gridSize = 10, progressInfo = TRUE){
  if(progressInfo) cat("Initialising profile CI\n")
  nParm <- length(coef(object))
  quant <- qchisq(p = level, df = 1) #df = length(coef(object)))
  llMod <- as.numeric(logLik(object))
  
  if(object$type == "binomial"){
    logLikNewParm <- function(par, object){
      n <- nrow(object$data)
      parmMat0 <- matrix(as.numeric(par), nrow = 1)
      
      prob <- object$pfFct(parmMat = parmMat0)(object$dataList$dose) 
      zeroTol <- 1e-8
      omZT <- 1 - zeroTol
      prob[prob > omZT] <- omZT
      prob[prob < zeroTol] <- zeroTol
      tmp <- -sum((object$dataList$resp * object$dataList$weights) * log(prob / (1 - prob)) + (object$dataList$weights * log(1 - prob)))
      
      total <- (object$"data")[, 5]
      success <- total*(object$"data")[, 2]    
      
      llVal <- sum(log(choose(total, success))) - tmp
      llVal
    }
  } 
  
  if(object$type == "continuous"){
    logLikNewParm <- function(par, object){
      n <- nrow(object$data)
      parmMat0 <- matrix(as.numeric(par), nrow = 1)
      
      fittedNewParm <- object$pfFct(parmMat = parmMat0)(object$dataList$dose)
      SSDNewParm <- sum((fittedNewParm - object$dataList$resp)^2)
      sigmaHatNewParm <- SSDNewParm/n
      
      llVal <- -n/2 * log(2*pi*sigmaHatNewParm) - SSDNewParm / (2*sigmaHatNewParm)
      llVal
    }
  }
  
  if(object$type == "count"){
    logLikNewParm <- function(par, object){
      n <- nrow(object$data)
      parmMat0 <- matrix(as.numeric(par), nrow = 1)
      
      fittedNewParm <- object$pfFct(parmMat = parmMat0)(object$dataList$dose)
      SSDNewParm <- sum((fittedNewParm - object$dataList$resp)^2)
      sigmaHatNewParm <- SSDNewParm/n
      
      llVal <- -n/2 * log(2*pi*sigmaHatNewParm) - SSDNewParm / (2*sigmaHatNewParm)
      llVal
    }
  }
  
  # Create grid
  confint0 <- confint(object, level = level)
  coef0 <- coef(object)
  
  gridPoints <- expand.grid(
    lapply(1:nParm, 
           function(i){ 
             seq(coef0[i] - 1.5 * (coef0[i] - confint0[i,1]), # Extra wide grids to ensure we cover entire confidence region
                 coef0[i] - 1.5 * (coef0[i] - confint0[i,2]), 
                 length.out = gridSize)
           }
    )
  )
  
  if(progressInfo){ 
    cat("Estimating confidence region\n")
    pb <- txtProgressBar(min = 0, max = nrow(gridPoints), style = 3)
  }
  llVals <- numeric(nrow(gridPoints))
  for(row in 1:nrow(gridPoints)){
    llVals[row] <- logLikNewParm(gridPoints[row,], object)
    if(progressInfo) setTxtProgressBar(pb, row)
  }
  if(progressInfo) close(pb)
  
  # Which values lie in the confidence region?
  accept <- 2 * (llMod - llVals) <= quant
  gridPointsAccept <- gridPoints[accept,]
  
  # Function for estimated bmd in new parameters
  bmdNewPar <- function(par){
    object0 <- object
    object0$fit$par <- par
    coefNames <- names(object0$coefficients)
    object0$coefficients <- par
    names(object0$coefficients) <- coefNames
    object0$curve[[1]] <- object$pfFct(parmMat = matrix(par, nrow = 1))
    object0$parmMat <- matrix(par, nrow = length(par))
    colnames(object0$parmMat) <- 1
    
    bmd(object0, bmr = bmr, backgType = backgType, backg = backg, controlSD = controlSD, def = def, interval = "delta", display = FALSE)
  }
  
  # Compute bmd values on all parameters in confidence region
  if(progressInfo){
    cat("Computing bmd values on parameters in confidence region\n")
    pb <- txtProgressBar(min = 0, max = sum(accept), style = 3)
  }
  bmdVals <- numeric(sum(accept))
  
  for(row in 1:sum(accept)){
    bmdVals[row] <- as.numeric(try(bmdNewPar(as.numeric(gridPointsAccept[row,]))$Results[1],
                                   silent = TRUE))
    if(progressInfo) setTxtProgressBar(pb, row)
  }
  if(progressInfo) close(pb)
  
  CI <- c(BMDL = min(bmdVals, na.rm = TRUE), BMDU = max(bmdVals, na.rm = TRUE))
  CI
}
