getCurveRepar <- function(x, backgType, background, controlSD, def, slope){
  model <- x
  
  if(!is.na(model$fct$fixed[4])){
    cat("Reparametrisation not defined for fixed parameter e.\n")
    return(NULL)
  }
  # Log-Logistic model
  else if (identical(substr(model$fct$name,1,2), "LL")){
    LL4repar <- function(BMD, par, bmr){
      # Handling fixed parameters first
      parVec <- model$fct$fixed
      parVec[is.na(parVec)] <- par
      b <- parVec[1]
      c <- parVec[2]
      d <- parVec[3]
      e <- parVec[4]
      
      # Background level
      if(identical(backgType, "modelBased")){
        p0 <- ifelse(identical(slope, "increasing"), c, d)
      } else {
        p0 <- background
      }
      
      if(identical(def, "excess")){
        z0 <- ifelse(identical(slope, "increasing"), (1-p0)*bmr + p0, -(1-p0)*bmr + p0)
      } else if(identical(def, "additional")){
        z0 <- ifelse(identical(slope, "increasing"), bmr + p0, -bmr + p0)
      } else if(identical(def, "point")){
        z0 <- bmr
      } else if(identical(def, "relative")){
        z0 <- ifelse(identical(slope, "increasing"), p0 + p0*bmr, p0 - p0*bmr)
      } else if(identical(def, "extra")){
        z0 <- ifelse(identical(slope, "increasing"), (d-p0)*bmr + p0, (c-p0)*bmr + p0)
      }
      e0 <- BMD * ( ( z0 - c) / (d - z0) )^(1/b)
      
      function(x){
        c + (d-c) / (1 + (x/e0)^b)
      }
    }
    LL4repar
  } 
  # Log-Normal model
  else if (identical(substr(model$fct$name,1,2), "LN")){
    LN4repar <- function(BMD, par, bmr){
      # Handling fixed parameters first
      parVec <- model$fct$fixed
      parVec[is.na(parVec)] <- par
      b <- parVec[1]
      c <- parVec[2]
      d <- parVec[3]
      e <- parVec[4]
      
      # Background level
      if(identical(backgType, "modelBased")){
        p0 <- ifelse(identical(slope, "increasing"), c, d)
      } else {
        p0 <- background
      }
      
      if(identical(def, "excess")){
        z0 <- ifelse(identical(slope, "increasing"), (1-p0)*bmr + p0, -(1-p0)*bmr + p0)
      } else if(identical(def, "additional")){
        z0 <- ifelse(identical(slope, "increasing"), bmr + p0, -bmr + p0)
      } else if(identical(def, "point")){
        z0 <- bmr
      } else if(identical(def, "relative")){
        z0 <- ifelse(identical(slope, "increasing"), p0 + p0*bmr, p0 - p0*bmr)
      } else if(identical(def, "extra")){
        z0 <- ifelse(identical(slope, "increasing"), (d-p0)*bmr + p0, (c-p0)*bmr + p0)
      }
      e0 <- BMD / exp( qnorm( (z0 - c) / (d-c) )/b)
      
      function(x){
        c + (d-c) * pnorm(b*(log(x) - log(e0)))
      }
    }
    LN4repar
  } 
  # Weibull 1 model
  else if (identical(substr(model$fct$name,1,2), "W1")){
    W14repar <- function(BMD, par, bmr){
      # Handling fixed parameters first
      parVec <- model$fct$fixed
      parVec[is.na(parVec)] <- par
      b <- parVec[1]
      c <- parVec[2]
      d <- parVec[3]
      e <- parVec[4]
      
      # Background level
      if(identical(backgType, "modelBased")){
        p0 <- ifelse(identical(slope, "increasing"), c, d)
      } else {
        p0 <- background
      }
      
      if(identical(def, "excess")){
        z0 <- ifelse(identical(slope, "increasing"), (1-p0)*bmr + p0, -(1-p0)*bmr + p0)
      } else if(identical(def, "additional")){
        z0 <- ifelse(identical(slope, "increasing"), bmr + p0, -bmr + p0)
      } else if(identical(def, "point")){
        z0 <- bmr
      } else if(identical(def, "relative")){
        z0 <- ifelse(identical(slope, "increasing"), p0 + p0*bmr, p0 - p0*bmr)
      } else if(identical(def, "extra")){
        z0 <- ifelse(identical(slope, "increasing"), (d-p0)*bmr + p0, (c-p0)*bmr + p0)
      }
      e0 <- BMD / (-log( (z0 - c) / (d - c)))^(1/b)
      
      function(x){
        c + (d-c) * exp(-exp(b*(log(x)-log(e0))))
      }
    }
    W14repar
  } 
  # Weibull 2 model
  else if (identical(substr(model$fct$name,1,2), "W2")){
    W24repar <- function(BMD, par, bmr){
      # Handling fixed parameters first
      parVec <- model$fct$fixed
      parVec[is.na(parVec)] <- par
      b <- parVec[1]
      c <- parVec[2]
      d <- parVec[3]
      e <- parVec[4]
      
      # Background level
      if(identical(backgType, "modelBased")){
        p0 <- ifelse(identical(slope, "increasing"), c, d)
      } else {
        p0 <- background
      }
      
      if(identical(def, "excess")){
        z0 <- ifelse(identical(slope, "increasing"), (1-p0)*bmr + p0, -(1-p0)*bmr + p0)
      } else if(identical(def, "additional")){
        z0 <- ifelse(identical(slope, "increasing"), bmr + p0, -bmr + p0)
      } else if(identical(def, "point")){
        z0 <- bmr
      } else if(identical(def, "relative")){
        z0 <- ifelse(identical(slope, "increasing"), p0 + p0*bmr, p0 - p0*bmr)
      } else if(identical(def, "extra")){
        z0 <- ifelse(identical(slope, "increasing"), (d-p0)*bmr + p0, (c-p0)*bmr + p0)
      }
      e0 <- BMD / (-log(1 - (z0 - c) / (d - c) ))^(1/b)
      
      
      function(x){
        c + (d-c) * (1-exp(-exp(b*(log(x)-log(e0)))))
      }
    }
    W24repar
  } 
  # Remaining models not reparametrised
  else{cat("Reparametrised curve not defined for model of type", model$fct$name, "\n")}
}

