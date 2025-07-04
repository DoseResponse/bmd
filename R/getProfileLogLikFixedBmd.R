getProfileLogLikFixedBmd <- function(object, curveRepar, bmr, start){
  n <- object$sumList$lenData
  dose <- object$dataList$dose
  response <- object$dataList$resp
  weights <- object$dataList$weights
  
  use_constrOptim <- TRUE
  # constrOptim options to ensure c < d when profiling likelihood (particularly important for Weibull models)
  if(all(is.na(object$fct$fixed[1:3]))){ # b,c,d not fixed
    ui <- matrix(c(0,-1,1), nrow = 1, byrow = TRUE)
    ci <- 0
  } else if(!is.na(object$fct$fixed[1]) & is.na(object$fct$fixed[2]) & is.na(object$fct$fixed[3])){ # b fixed
    ui <- matrix(c(-1,1), nrow = 1)
    ci <- 0 
  } else if(!is.na(object$fct$fixed[1]) & !is.na(object$fct$fixed[2]) & is.na(object$fct$fixed[3])){ # b,c fixed
    ui <- matrix(1, nrow = 1)
    ci <- object$fct$fixed[2]
  } else if(!is.na(object$fct$fixed[1]) & is.na(object$fct$fixed[2]) & !is.na(object$fct$fixed[3])){ # b,d fixed
    ui <- matrix(-1, nrow = 1)
    ci <- -object$fct$fixed[3]
  } else if(is.na(object$fct$fixed[1]) & !is.na(object$fct$fixed[2]) & is.na(object$fct$fixed[3])){ # c fixed
    ui <- matrix(c(0,1), nrow = 1)
    ci <- object$fct$fixed[2]
  } else if(is.na(object$fct$fixed[1]) & !is.na(object$fct$fixed[2]) & !is.na(object$fct$fixed[3])){ # c,d fixed
    use_constrOptim <- FALSE
    ui <- matrix(0, nrow = 1)
    ci <- 0
  } else if(is.na(object$fct$fixed[1]) & is.na(object$fct$fixed[2]) & !is.na(object$fct$fixed[3])){ # d fixed
    ui <- matrix(c(0,-1), nrow = 1)
    ci <- -object$fct$fixed[3]
  }
  
  # define profile log-likelihood
  if(identical(object$type, "continuous")){
    profileLogLikFixedBmd <- function(BMD){
      fn0 <- function(par){sum(((response - curveRepar(BMD, par, bmr)(dose))/weights)^2)}
      if(use_constrOptim){
        constrOptim0 <- constrOptim(theta=start, 
                                  f = fn0,
                                  grad = NULL,
                                  ui = ui,
                                  ci = ci)
        SSD <- constrOptim0$value
      } else {
        int0 <- confint(object,level = 0.99)[1,]
        optim0 <- optim(par = start, fn = fn0, method = "Brent", lower = int0[1], upper = int0[2])
        SSD <- optim0$value
      }
      #sigmaSqHat <- object$fit$value/n  # 
      
      #llVal <- -n/2 * log(2*pi*sigmaSqHat) - SSD / (2*sigmaSqHat)
      llVal <- -(n/2)*(log(2*pi)+log(SSD)-log(n)+1)
      llVal
    }
  } else if(identical(object$type, "binomial")){
    profileLogLikFixedBmd <- function(BMD){
      zeroTol <- 1e-8
      fn0 <- function(par)  # dose, resp and weights are fixed
      {
        prob <- curveRepar(BMD, par, bmr)(dose) # multCurves(dose / doseScaling, c)
        omZT <- 1 - zeroTol
        prob[prob > omZT] <- omZT
        prob[prob < zeroTol] <- zeroTol
        -sum((response * weights) * log(prob / (1 - prob)) + (weights * log(1 - prob)))
      }
      
      constrOptim0 <- constrOptim(theta=start, 
                                  f = fn0,
                                  grad = NULL,
                                  ui = ui,
                                  ci = ci)
      
      total <- (object$"data")[, 5]
      success <- total*(object$"data")[, 2]    
      
      llVal <- sum(log(choose(total, success))) - constrOptim0$value
      llVal
    }
  } else if(identical(object$type, "Poisson")){
    profileLogLikFixedBmd <- function(BMD){
      fn0 <- function(par){sum(response * log(curveRepar(BMD, par, bmr)(dose)) - curveRepar(BMD, par, bmr)(dose) )} # sum(response * log(curveRepar(BMD, par, bmr)(dose)) - curveRepar(BMD, par, bmr)(dose) - log(foctorial(response)) )
      constrOptim0 <- constrOptim(theta=start, 
                                  f = fn0,
                                  grad = NULL,
                                  ui = ui,
                                  ci = ci)
      
      llVal <- constrOptim0$value
      llVal
    }
  }
  
  profileLogLikFixedBmd
}
