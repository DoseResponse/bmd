# ED.bmd <- function(object, 
#                    respLev, 
#                    interval = c("none", "delta", "fls", "tfls", "inv"), 
#                    clevel = NULL,
#                    level = ifelse(!(interval == "none"), 0.95, NULL), 
#                    reference = c("control", "upper"), 
#                    type = c("relative", "absolute"), 
#                    lref, uref, bound = TRUE, 
#                    vcov. = vcov,
#                    display = TRUE, 
#                    logBase = NULL, 
#                    multcomp = FALSE, 
#                    intType = "confidence", ...)
# {
#   interval <- match.arg(interval)
#   reference <- match.arg(reference)
#   type <- match.arg(type)
#   
#   ## Checking 'respLev' vector: it should be numbers between 0 and 100
#   if ( (type == "relative") && (bound) ) 
#   {
#     if (any(respLev <= 0 | respLev >= 100)) 
#     {
#       stop("Response levels (percentages) outside the interval ]0, 100[ not allowed")
#     }
#   }
#   
#   ## Retrieving relevant quantities
#   EDlist <- object$fct[["edfct"]]  
#   if (is.null(EDlist)) {stop("ED values cannot be calculated")}         
#   indexMat <- object[["indexMat"]]
#   parmMat <- object[["parmMat"]]
#   
#   curveNames <- colnames(parmMat)  # colnames(object$"parmMat")
#   options(warn = -1)  # switching off warnings caused by coercion in the if statement
#   if (any(is.na(as.numeric(curveNames))))
#   {
#     curveOrder <- order(curveNames)
#   } else { # if names are numbers then skip re-ordering
#     curveOrder <- 1:length(curveNames)
#   }
#   options(warn = 0)  # normalizing behaviour of warnings
#   
#   strParm0 <- curveNames[curveOrder]
#   indexMat <- indexMat[, curveOrder, drop = FALSE]
#   parmMat <- parmMat[, curveOrder, drop = FALSE]
#   
#   strParm <- strParm0
#   #vcMat <- vcov.(object)
#   if (is.function(vcov.))  # following a suggestion by Andrea Onofri
#   {
#     vcMat <- vcov.(object)
#   } else {
#     vcMat <- vcov.
#   }
#   
#   ## Defining vectors and matrices needed
#   ncolIM <- ncol(indexMat)
#   indexVec <- 1:ncolIM    
#   #    lenEB <- ncolIM    
#   lenPV <- length(respLev)  # used twice below
#   noRows <- ncolIM * lenPV
#   dimNames <- rep("", noRows)  # lenEB*lenPV, 2)
#   EDmat <- matrix(0, noRows, 2)  # lenEB*lenPV, 2)
#   oriMat <- matrix(0, noRows, 2)  # lenEB*lenPV, 2)
#   pmodelsMatrixList <- getpmodelsMatrixList(object)
#   
#   ## Skipping curve id if only one curve is present
#   if (identical(length(unique(strParm)), 1)) 
#   {
#     strParm[indexVec] <- rep("", ncolIM)
#   } else {
#     strParm <- paste(strParm, ":", sep = "")
#   }
#   
#   ## Calculating estimates and estimated standard errors
#   rowIndex <- 1
#   lenIV <- length(indexVec)
#   # dEDmat <- matrix(0, lenPV * lenIV, nrow(vcMat))
#   intMat <- NULL
#   for (i in indexVec)
#   {
#     parmChosen <- parmMat[, i]
#     #parmInd <- indexMat[, i]
#     #varCov <- vcMat[parmInd, parmInd]
#     varCov <- pmodelsMatrixList[[curveNames[i]]] %*% vcMat %*% t(pmodelsMatrixList[[curveNames[i]]])
#     
#     if ((is.null(clevel)) || (strParm0[i] %in% clevel))
#     {
#       for (j in 1:lenPV)
#       {
#         EDeval <- EDlist(parmChosen, respLev[j], reference = reference, type = type, ...)            
#         EDval <- EDeval[[1]]
#         dEDval <- EDeval[[2]]
#         # dEDmat[(i-1)*lenPV + j, parmInd] <- dEDval 
#         
#         oriMat[rowIndex, 1] <- EDval
#         oriMat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
#         
#         if (!is.null(logBase))
#         {
#           EDval <- logBase^(EDval)                
#           dEDval <- EDval * log(logBase) * dEDval
#         }
#         EDmat[rowIndex, 1] <- EDval
#         EDmat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
#         
#         dimNames[rowIndex] <- paste(strParm[i], respLev[j], sep = "")
#         rowIndex <- rowIndex + 1
#       }
#       if (interval == "inv")
#       {
#         intMat <- rbind(intMat, t(EDinvreg1(object, respLev, strParm0[i], 
#                                             intType = intType, level = level, type = type)))
#       }
#       
#     } else {
#       rowsToRemove <- rowIndex:(rowIndex + lenPV - 1)
#       EDmat <- EDmat[-rowsToRemove, , drop = FALSE]
#       dimNames <- dimNames[-rowsToRemove]
#     }
#     
#   }
#   
#   ## Defining column names
#   colNames <- c("Estimate", "Std. Error")
#   
#   ## Calculating the confidence intervals
#   if (interval == "delta")
#   {
#     intMat <- drc:::confint.basic(EDmat, level, object$"type", df.residual(object), FALSE)
#     intLabel <- "Delta method"
#   }
#   
#   if (interval == "tfls")
#   {
#     intMat <- exp(drc:::confint.basic(matrix(c(log(oriMat[, 1]), oriMat[, 2] / oriMat[, 1]), ncol = 2), 
#                                       level, object$"type", df.residual(object), FALSE))
#     intLabel <- "To and from log scale"       
#   }
#   
#   if (interval == "fls")
#   {        
#     if (is.null(logBase)) 
#     {
#       logBase <- exp(1)
#       EDmat[, 1] <- exp(EDmat[, 1])  # back-transforming log ED values
#     }
#     
#     intMat <- logBase^(drc:::confint.basic(oriMat, level, object$"type", df.residual(object), FALSE))
#     intLabel <- "Back-transformed from log scale"  
#     
#     ## Dropping estimated standard errors (not relevant after back transformation)        
#     EDmat <- EDmat[, -2, drop = FALSE]  
#     colNames <- colNames[-2]
#     #        colNames <- c(colNames[-2], "Lower", "Upper")  # standard errors not relevant        
#   }
#   
#   if (interval == "inv")
#   {
#     EDmat <- EDmat[, -2, drop = FALSE]  
#     colNames <- colNames[-2]
#     intLabel <- "Inverse regression"  
#   }
#   
#   if (identical(interval, "none"))
#   {
#     intLabel <- NULL
#   } else {
#     EDmat <- as.matrix(cbind(EDmat, intMat))
#     colNames <- c(colNames, "Lower", "Upper")         
#   } 
#   dimnames(EDmat) <- list(dimNames, colNames)
#   rownames(EDmat) <- paste("e", rownames(EDmat), sep = ":")
#   
#   EDmat      
# }
