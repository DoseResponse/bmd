bmd.edfct <- function(object){
  if(class(object$fct) %in% c("llogistic", "log-normal", "Weibull-1", "Weibull-2", "Boltzmann", "braincousens", "fp-logistic")){
    ## Handling 'fixed' argument
    numParm <- length(object$fct$fixed)
    notFixed <- is.na(object$fct$fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- object$fct$fixed[!notFixed]
    
    # Log-logistic
    if(identical(class(object$fct), "llogistic")){
      if(substr(object$fct$name, 3,3) == "."){
        edfct <- function(parm, respl, reference, type, ...)
        {
          parmVec[notFixed] <- parm
          p <- EDhelper(parmVec, respl, reference, type)
          
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
      } else if(substr(object$fct$name, 3,3) == "2"){
        edfct <- function(parm, respl, reference, type, ...)
        {
          parmVec[notFixed] <- parm
          
          p <- EDhelper(parmVec, respl, reference, type)
          
          tempVal <- log((100-p)/100)
          EDp <- exp(parmVec[4])*(exp(-tempVal/parmVec[5])-1)^(1/parmVec[1])
          
          EDder <- 
            EDp*c(-log(exp(-tempVal/parmVec[5])-1)/(parmVec[1]^2), 
                  0, 0, 1, 
                  exp(-tempVal/parmVec[5])*tempVal/(parmVec[5]^2)*(1/parmVec[1])*((exp(-tempVal/parmVec[5])-1)^(-1)))
          
          if (type == "absolute")
          {
            EDder[2:3] <- EDp * exp(-tempVal/parmVec[5]) / (parmVec[1]*parmVec[5]*(parmVec[3]-parmVec[2])*(exp(-tempVal/parmVec[5])-1)) *
              c( exp(-tempVal)-1, 1)
          }
          
          return(list(EDp, EDder[notFixed]))
        }
      }
    }
    
    # Log-Normal
    if(identical(class(object$fct), "log-normal")){
      edfct <- function(parm, respl, reference, type, ...)
      {
        parmVec[notFixed] <- parm
        p <- absToRel(parmVec, respl, type)
        
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
        p <- EDhelper(parmVec, respl, reference, type)
        
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
        
        p <- absToRel(parmVec, respl, type)
        
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
    
    if(identical(class(object$fct), "Boltzmann")){
      edfct <- function(parm, respl, reference, type, ...)
      {
        parmVec[notFixed] <- parm
        p <- EDhelper(parmVec, respl, reference, type)
        
        #        if (parmVec[1] > 0) 
        #        {
        #            tempVal <- (100 - p) / 100
        # old=wrong            EDp <- parmVec[4] + log( (1 + exp(-parmVec[1]*parmVec[4])) / (tempVal^(1/parmVec[5])) - 1)/parmVec[1]
        
        # old=wrong            ## deriv(~e + log( (1 + exp(-b*e)) / (((100 - p) / 100)^(1/f)) - 1)/b, c("b", "c", "d", "e", "f"), function(b,c,d,e,f){})
        
        ## deriv(~e + log((100/(100-p))^(1/f) - 1) / b, c("b", "c", "d", "e", "f"), function(b,c,d,e,f){})
        ## evaluated at the R prompt
        EDderFct <- 
          function (b, c, d, e, f) 
          {
            .expr2 <- 100/(100-p)
            .expr4 <- .expr2^(1/f)
            .expr5 <- .expr4 - 1
            .expr6 <- log(.expr5)
            .value <- e + .expr6/b
            .grad <- array(0, c(length(.value), 5L), list(NULL, c("b", "c", "d", "e", "f")))
            .grad[, "b"] <- -(.expr6/b^2)
            .grad[, "c"] <- ifelse(identical(type, "absolute"),
                                   .expr4 * (.expr2-1) / (b*f*(d-c)*.expr5),
                                   0)
            .grad[, "d"] <- ifelse(identical(type, "absolute"), 
                                   .expr4 / (b*f*(d-c)*.expr5),
                                   0)
            .grad[, "e"] <- 1
            .grad[, "f"] <- -(.expr4 * (log(.expr2) * (1/f^2))/.expr5/b)
            attr(.value, "gradient") <- .grad
            .value
          }
        EDcalc <- EDderFct(parmVec[1], parmVec[2], parmVec[3], parmVec[4], parmVec[5])
        EDp <- as.numeric(EDcalc)
        EDder <- attr(EDcalc, "gradient")
    
        return(list(EDp, EDder[notFixed]))
      }
    }
    
    if(identical(class(object$fct), "braincousens")){
      edfct <- function(parm, respl, reference, type, lower = 1e-3, upper = 10000, ...)
      {
        #        if (is.missing(upper)) {upper <- 1000}
        interval <- c(lower, upper)     
        
        parmVec[notFixed] <- parm
        p <- EDhelper(parmVec, respl, reference, type)
        tempVal <- (100-p)/100
        
        edfct0 <- function(parmVec){
          p <- EDhelper(parmVec, respl, reference, type)
          tempVal <- (100-p)/100
          
          helpEqn <- function(dose) 
          {
            expVal <- exp(parmVec[1]*(log(dose)-log(parmVec[4])))
            parmVec[5]*(1+expVal*(1-parmVec[1]))-(parmVec[3]-parmVec[2])*expVal*parmVec[1]/dose
          }
          maxAt <- uniroot(helpEqn, interval)$root
          
          eqn <- function(dose) {tempVal*(1+exp(parmVec[1]*(log(dose)-log(parmVec[4]))))-(1+parmVec[5]*dose/(parmVec[3]-parmVec[2]))}
          EDp <- uniroot(eqn, lower = maxAt, upper = upper)$root
          EDp
        }
        
        EDp <- edfct0(parmVec)
        EDdose <- EDp
        tempVal1 <- exp(parmVec[1]*(log(EDdose)-log(parmVec[4])))
        tempVal2 <- parmVec[3]-parmVec[2]
        derParm <- c(tempVal*tempVal1*(log(EDdose)-log(parmVec[4])), -parmVec[5]*EDdose/((tempVal2)^2),
                     parmVec[5]*EDdose/((tempVal2)^2), -tempVal*tempVal1*parmVec[1]/parmVec[4],
                     -EDdose/tempVal2)
        derParm <- c(-tempVal*tempVal1*(log(EDdose)-log(parmVec[4])), parmVec[5]*EDdose/((tempVal2)^2),
                     parmVec[5]*EDdose/((tempVal2)^2), tempVal*tempVal1*parmVec[1]/parmVec[4],
                     EDdose/tempVal2)
        derDose <- tempVal*tempVal1*parmVec[1]/EDdose-parmVec[5]/tempVal2
        
        EDder <- derParm/derDose
        EDder[2:3] <- numDeriv::jacobian(edfct0, parmVec)[2:3]
        
        return(list(EDp, EDder[notFixed]))
      }
    }
    
    if(identical(class(object$fct), "fp-logistic")){
      if(!requireNamespace("numDeriv")){
        stop('package "numDeriv" must be installed to use FPL models')
      }
      
      p1 <- as.numeric(unlist(strsplit(object$fct$name, split = "[,()]+"))[2])
      p2 <- as.numeric(unlist(strsplit(object$fct$name, split = "[,()]+"))[3])
      edfct <- function(parm, respl, reference, type, loged = FALSE, ...)
      {
        parmVec[notFixed] <- parm
        p <- EDhelper2(parmVec, respl, reference, type, parmVec[1] > 0)
        
        invfp <- function(resp, b, e)
        {
          fct0 <- function(x){resp - (b*(log(x+1)^p1) + e*(log(x+1)^p2))}
          uniroot(fct0, c(0.001, 5000))$root
        }
        
        EDfct <- function(b, c, d, e) 
        {
          invfp(log((100-p)/p), b, e)
        }
        
        EDfct0 <- function(par) 
        {
          p <- EDhelper2(par, respl, reference, type, par[1] > 0)
          invfp(log((100-p)/p), par[1], par[4])
        }
        
        EDp <- EDfct(parmVec[1], parmVec[2], parmVec[3], parmVec[4])
        
        logEDp <- log(EDp+1)
        denVal <- parmVec[1] * p1 * (logEDp)^(p1-1) + parmVec[4] * p2 * (logEDp)^(p2-1)
        derVec <- (EDp+1) * c(logEDp^p1, logEDp^p2) / denVal
        EDder <- -c(derVec[1], 0, 0, derVec[2])
        
        # Addition by Jens Riis Baalkilde
        if(type == "absolute"){
          EDder[2:3] <- numDeriv::jacobian(EDfct0, parmVec)[2:3]
        }
        
        # if(loged) 
        # {
        #   EDder <- EDder / EDp
        #   EDp <- log(EDp)
        # }
        return(list(EDp, EDder[notFixed]))
      }
    }
  } else {
    edfct <- object$fct[["edfct"]]
  }
  
  edfct
}
