# HELPER FUNCTION BORROWED FROM DRC PACKAGE
# Defining basic function for providing confidence intervals

"confint.basic" <- function(estMat, level, intType, dfres, formatting = TRUE)
{
  alphah <- (1 - level)/2 
  #    if (type == "u") {two <- qnorm(1 - alphah)}
  #    if (type == "t") {two <- qt(1 - alphah, df.residual(object))}    
  tailPercentile <- switch(intType, 
                           binomial = qnorm(1 - alphah), 
                           continuous = qt(1 - alphah, dfres),
                           event = qnorm(1 - alphah),
                           Poisson = qnorm(1 - alphah),
                           negbin1 = qnorm(1 - alphah),
                           negbin2 = qnorm(1 - alphah))
  
  estVec <- estMat[, 1]
  halfLength <- tailPercentile * estMat[, 2]
  confMat <- matrix(c(estVec -  halfLength, estVec +  halfLength), ncol = 2)    
  
  ## Formatting matrix
  if (formatting)
  {
    colnames(confMat) <- c(paste(format(100 * alphah), "%", sep = " "), paste(format(100*(1 - alphah)), "%", sep = " "))
    rownames(confMat) <- rownames(estMat)
  }
  
  return(confMat)    
}



