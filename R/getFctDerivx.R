getFctDerivx <- function(object){
  ## Handling 'fixed' argument
  numParm <- length(object$fct$fixed)
  notFixed <- is.na(object$fct$fixed)
  parmVec <- rep(0, numParm)
  parmVec[!notFixed] <- object$fct$fixed[!notFixed]
  
  derivx <- NULL
  # Log-logistic
  if(identical(class(object$fct), "llogistic")){
    if(substr(object$fct$name, 3,3) == "."){
      derivx <- function(x, parm)
      {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        #        bNeg <- parmMat[, 1] < 0
        #        parmMat[bNeg, 1] <- -parmMat[bNeg, 1]
        
        temp1 <- x/parmMat[, 4]          
        temp2 <- 1 + (temp1)^parmMat[, 1]
        temp3 <- parmMat[, 5]*(temp2^(parmMat[, 5] - 1))*(parmMat[, 1]/parmMat[, 4])*temp1^(parmMat[, 1] - 1)
        temp4 <- temp2^(2*parmMat[, 5])
        
        (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4
        retVec <- (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4
        #        retVec[bNeg] <- -retVec[bNeg]
        retVec
      }
    } else if(substr(object$fct$name, 3,3) == "2"){
      derivx <- function(x, parm)
      {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        
        temp1 <- exp(parmMat[, 1]*(log(x) - parmMat[, 4]))  # x/parmMat[, 4]          
        temp2 <- 1 + temp1
        temp3 <- parmMat[, 5]*(temp2^(parmMat[, 5] - 1))*temp1*parmMat[, 1]/x
        temp4 <- temp2^(2*parmMat[, 5])
        
        (-(parmMat[, 3] - parmMat[, 2])*temp3)/temp4 
      }
    }
  }
  
  # Log-Normal
  if(identical(class(object$fct), "log-normal")){
    derivx <- function(dose, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      
      dFct <- function (dose, b, c, d, e) 
      {
        .expr1 <- d - c
        .expr5 <- b * (log(dose) - transfe(e))
        .value <- c + .expr1 * pnorm(.expr5)
        .grad <- array(0, c(length(.value), 1L), list(NULL, c("dose")))
        .grad[, "dose"] <- .expr1 * (dnorm(.expr5) * (b * (1/dose)))
        attr(.value, "gradient") <- .grad
        .value
      }        
      attr(dFct(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4]), "gradient")
    }
  }
  
  # Weibull1
  if(identical(class(object$fct), "Weibull-1")){
    derivx <- function(x, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      
      .expr1 <- parmMat[, 3] - parmMat[, 2]  # d - c
      .expr6 <- exp(parmMat[, 1] * (log(x) - log(parmMat[, 4])))
      .expr8 <- exp(-.expr6)
      .value <- parmMat[, 2] + .expr1 * .expr8
      .grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
      .grad[, "x"] <- -(.expr1 * (.expr8 * (.expr6 * (parmMat[, 1] * (1/x)))))
      .grad
    }
  }
  
  # Weibull2
  if(identical(class(object$fct), "Weibull-2")){
    derivx <- function(x, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      
      .expr1 <- parmMat[, 3] - parmMat[, 2]
      .expr6 <- exp(parmMat[, 1] * (log(x) - log(parmMat[, 4])))
      .expr8 <- exp(-.expr6)
      .value <- parmMat[, 2] + .expr1 * (1 - .expr8)
      .grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
      .grad[, "x"] <- .expr1 * (.expr8 * (.expr6 * (parmMat[, 1] * (1/x))))
      .grad
    }
  }
  
  if(identical(class(object$fct), "Boltzmann")){
    derivx <- function(x, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      
      temp1 <- exp(parmMat[, 1]*(x - parmMat[, 4]))       
      
      (-parmMat[, 5]*(parmMat[, 3] - parmMat[, 2])*temp1*parmMat[, 1])/((1 + temp1)^(parmMat[, 5] + 1))
    }
  }
  
  # Brain-Cousens
  if(identical(class(object$fct), "braincousens")){
    derivx <- function(dose, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      
      t1 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
      t2 <- 1 + t1
      
      parmMat[, 5] / t2 - (parmMat[, 5]*dose - parmMat[, 2] + parmMat[, 3]) * t1 * parmMat[, 1] / (t2^2 * dose)
    }
  }
  
  
  if(identical(class(object$fct), "fp-logistic")){
    p1 <- as.numeric(unlist(strsplit(object$fct$name, split = "[,()]+"))[2])
    p2 <- as.numeric(unlist(strsplit(object$fct$name, split = "[,()]+"))[3])
    derivx <- function(dose, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      
      dFct <- function (dose, b, c, d, e) 
      {
        .expr1 <- d - c
        .expr2 <- dose + 1
        .expr3 <- log(.expr2)
        .expr9 <- exp(b * .expr3^p1 + e * .expr3^p2)
        .expr10 <- 1 + .expr9
        .expr15 <- 1/.expr2
        .value <- c + .expr1/.expr10
        .grad <- array(0, c(length(.value), 1L), list(NULL, c("dose")))
        .grad[, "dose"] <- -(.expr1 * (.expr9 * (b * (.expr3^(p1 - 1) * (p1 * .expr15)) + e * (.expr3^(p2 - 1) * (p2 * .expr15))))/.expr10^2)
        attr(.value, "gradient") <- .grad
        .value
      }
      attr(dFct(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4]), "gradient")
    }
  }

  derivx
}
