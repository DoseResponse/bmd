bmd.edfct <- function(object){
  if(class(object$fct) %in% c("llogistic", "log-normal", "Weibull-1", "Weibull-2")){
    ## Handling 'fixed' argument
    numParm <- length(object$fct$fixed)
    notFixed <- is.na(object$fct$fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- object$fct$fixed[!notFixed]
    
    # Log-logistic
    if(identical(class(object$fct), "llogistic")){
      edfct <- function(parm, respl, reference, type, ...)
      {
        parmVec[notFixed] <- parm
        p <- drc:::EDhelper(parmVec, respl, reference, type)
        
        tempVal <- log((100-p)/100)
        EDp <- parmVec[4]*(exp(-tempVal/parmVec[5])-1)^(1/parmVec[1])
        
        EDder <- 
          EDp*c(-log(exp(-tempVal/parmVec[5])-1)/(parmVec[1]^2), 
                0, 0, 1/parmVec[4], 
                exp(-tempVal/parmVec[5])*tempVal/(parmVec[5]^2)*(1/parmVec[1])*((exp(-tempVal/parmVec[5])-1)^(-1)))
        
        if (type == "absolute")
        {
          EDder[2:3] <- EDp * exp(-tempVal/parmVec[5]) / (parmVec[1]*parmVec[5]*(parmVec[3]-parmVec[2])*(exp(-tempVal/parmVec[5])-1)) *
            c( exp(-tempVal)-1, 1)
        }
        
        return(list(EDp, EDder[notFixed]))
      }
    }
    
    # Log-Normal
    if(identical(class(object$fct), "log-normal")){
      edfct <- function(parm, respl, reference, type, ...)
      {
        parmVec[notFixed] <- parm
        p <- drc:::absToRel(parmVec, respl, type)
        
        ## Reversing p
        if (identical(type, "absolute"))
        {
          p <- 100 - p
        }
        if (identical(type, "relative") && (parmVec[1] < 0) && (reference == "control"))
        {
          p <- 100 - p
        }
        
        pProp <- 1 - (100-p) / 100
        
        ## deriv(~e * exp(22 / b), c("b", "c", "d", "e"), function(b,c,d,e){})
        ## using "22" instead of qnorm(pProp)
        EDfct <- function (b, c, d, e) 
        {
          .expr2 <- exp(qnorm(pProp) / b)
          .value <- e * .expr2
          .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
          .grad[, "b"] <- -(e * (.expr2 * (qnorm(pProp) / (b^2))))
          .grad[, "c"] <- ifelse(identical(type, "absolute"), 
                                 .value * (pProp-1)/(b*dnorm(qnorm(pProp))*(d-c)),
                                 0)  # 0
          .grad[, "d"] <- ifelse(identical(type, "absolute"), 
                                 - .value * pProp/(b*dnorm(qnorm(pProp))*(d-c)),
                                 0) # 0
          .grad[, "e"] <- .expr2
          attr(.value, "gradient") <- .grad
          .value
        }
        
        EDp <- EDfct(parmVec[1], parmVec[2], parmVec[3], parmVec[4])
        EDder <- attr(EDfct(parmVec[1], parmVec[2], parmVec[3], parmVec[4]), "gradient")
        return(list(EDp, EDder[notFixed]))
      }
    }
    
    # Weibull1
    if(identical(class(object$fct), "Weibull-1")){
      edfct <- function(parm, respl, reference, type, ...)  # function(parm, p, reference, type, ...)
      {        
        parmVec[notFixed] <- parm
        p <- drc:::EDhelper(parmVec, respl, reference, type)
        
        tempVal <- log(-log((100-p)/100))
        EDp <- exp(tempVal/parmVec[1] + log(parmVec[4]))
        
        EDder <- EDp*c(-tempVal/(parmVec[1]^2), 0, 0, 1/parmVec[4])
        
        if(type == "absolute"){
          tempVal2 <- (100-p)/100
          EDder[2:3] <- c(EDp * (tempVal2 - 1)/(parmVec[1] * tempVal2 * (parmVec[3]-parmVec[2]) * log(tempVal2)),
                          - EDp / (parmVec[1] * (parmVec[3]-parmVec[2]) * log(tempVal2)) )
        }
        
        return(list(EDp, EDder[notFixed]))
      }
    }
    
    # Weibull2
    if(identical(class(object$fct), "Weibull-2")){
      edfct <- function(parm, respl, reference, type, ...)
      {   
        parmVec[notFixed] <- parm
        
        p <- drc:::absToRel(parmVec, respl, type)
        
        ## Reversing p
        if ( (parmVec[1] > 0) && (reference == "control") && (type == "relative") )
        {
          p <- 100 - p
        }
        
        tempVal <- log(-log(p/100))
        EDp <- exp(tempVal/parmVec[1] + log(parmVec[4]))
        EDder <- EDp*c(-tempVal/(parmVec[1]^2), 0, 0, 1/parmVec[4])
        
        if(identical(type, "absolute")){
          tempVal2 <- p/100
          EDder[2:3] <- c(EDp * tempVal2/(parm[1] * tempVal2 * (parm[3]-parm[2]) * log(tempVal2)),
                          EDp * (1-tempVal2)/ (parm[1] * (parm[3]-parm[2]) * tempVal2 * log(tempVal2)) )
        }
        
        return(list(EDp, EDder[notFixed]))
      }
    }
  } else {
    edfct <- object$fct[["edfct"]]
  }
  
  edfct
}
