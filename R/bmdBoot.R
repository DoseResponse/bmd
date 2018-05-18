bmdBoot <- function(object, bmr, R=1000, boot="resample", bmdType = "orig",
                    backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
                    backg=NA, 
                    def = c("excess", "additional", 
                            "relative", "extra", "added", "hybridExc", "hybridAdd", "point")){
  if(boot=="resample"){
    if(object$type=="binomial"){
      data.str <- object$data
      data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
      data.e<-expandBinomial(data.str, 
                          number = "number",
                          total = "weights",
                          dose = as.character(object$call$formula[[3]]))
  data.e[,"row.num"]<-1:dim(data.e)[1]
  tmp.data <- list()
  for(i in 1:R){
    sampled.expand <- data.e[as.numeric(unlist(aggregate(row.num ~ data.e[,as.character(object$call$formula[[3]])], data=data.e, 
                                                        FUN=function(x) sample(x,replace=TRUE))[[2]])),]
    df <- aggregate(cbind(sampled.expand[,"number"],
                          sampled.expand[,"weights"]) ~ 
                      sampled.expand[,as.character(object$call$formula[[3]])],FUN = sum)
    colnames(df) <- c(as.character(object$call$formula[[3]]),
                      as.character(object$call$formula[[2]])[[2]],
                      as.character(object$call$formula[[2]])[[3]])
    tmp.data[[i]] <- df
  }
    }
    if(object$type=="continuous"){
      data.e<-object$data
      data.e[,"row.num"]<-1:dim(data.e)[1]
      data.e[,"dose"]<-data.e[,as.character(object$call$formula[[3]])]
      tmp.data <- list()
      for(i in 1:R){
        tmp.data[[i]] <- data.e[as.numeric(unlist(aggregate(row.num ~ dose, data=data.e, 
                                                             FUN=function(x) sample(x,replace=TRUE))[[2]])),]
         }
      }
  } else if(boot=="pseudorandom"){
    if(object$type=="binomial"){
    Y <- object$data[[as.character(object$call$formula)[[2]]]]*object$data[["weights"]]
    N <- object$data[["weights"]]
    shrinks <- which(Y==N | Y==0)
    Y[shrinks] <- Y[shrinks]+0.25
    N[shrinks] <- N[shrinks]+0.5
    prob <- rep(Y/N,N)
    tmp.data <- list()
    for(i in 1:R){
      sampled.expand <- data.frame(number = rbinom(length(prob),1,prob), 
                                   dose = rep(object$data[,as.character(object$call$formula[[3]])],N), 
                                   total = 1)
      df <- aggregateBinomial(number/total~dose, sampled.expand)
      colnames(df) <- c(as.character(object$call$formula[[3]]),
                        as.character(object$call$formula[[2]])[[2]],
                        as.character(object$call$formula[[2]])[[3]])
      tmp.data[[i]] <- df  
    }
    }
    if(object$type=="continuous"){
      mean.Y <- aggregate(object$data[,as.character(object$call$formula[[2]])] ~ 
                            object$data[,as.character(object$call$formula[[3]])],
                          FUN=function(x) mean(x,na.rm=TRUE))[,2]
      sd.Y <- aggregate(object$data[,as.character(object$call$formula[[2]])] ~ 
                          object$data[,as.character(object$call$formula[[3]])],
                        FUN=function(x) sd(x,na.rm=TRUE))[,2]
      Dose<- aggregate(object$data[,as.character(object$call$formula[[2]])] ~ 
                         object$data[,as.character(object$call$formula[[3]])],
                       FUN=function(x) length(!is.na(x)))[,1]
      N.dose<- aggregate(object$data[,as.character(object$call$formula[[2]])] ~ 
                           object$data[,as.character(object$call$formula[[3]])],
                         FUN=function(x) length(!is.na(x)))[,2]
      tmp.data <- list()
      for(i in 1:R){
        sampled <- data.frame(y = rnorm(sum(N.dose),mean=rep(mean.Y,N.dose),sd=rep(sd.Y,N.dose)), 
                                     dose = rep(Dose,N.dose))
        colnames(sampled) <- c(as.character(object$call$formula[[2]]), as.character(object$call$formula[[3]]))
        tmp.data[[i]] <- sampled
        }
    }
  }
    drm.list <- lapply(tmp.data, function(x){
      drm(object$call$formula, data = x, type = object$type, fct = object[["fct"]])}
      )
      
    bmd.list <- lapply(drm.list,function(x){
    bmd(x, bmr = bmr, backgType = backgType, backg=backg, def=def)[1]}
    )
    if(bmdType == "orig"){
      use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def)[1]
    } else if(bmdType == "mean"){
      use.bmd <- mean(unlist(bmd.list))
    } else if(bmdType == "median"){
      use.bmd <- quantile(unlist(bmd.list),c(0.5))  
      } 
    resMat <- matrix(NA,1,2)
    resMat[1,1] <- use.bmd
    resMat[1,2] <- quantile(unlist(bmd.list),c(0.05))
    colnames(resMat) <- c("BMD", "BMDL")
    rownames(resMat) <- c("")
    resMat    
   
}


