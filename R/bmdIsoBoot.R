bmdIsoBoot <- function(object, data, type, bmr, R=1000, boot="resample", 
                       backgType = c("modelBased", "absolute","hybridSD","hybridPercentile"), 
                       backg=NA, 
                       def = c("excess", "additional", "relative", "added", "hybridExc", "hybridAdd", "point")){
  if (type %in% c("Poisson","negbin1","negbin2") & boot!="resample") {
    stop(paste("\"",type,"\" only works with resample bootstrap \"", sep=""))
  }
  if(boot=="resample"){
    if(type=="binomial"){
      data.e<-expandBinomial(data, 
                          number = strsplit(as.character(object[[2]]),"/")[[2]],
                          total = strsplit(as.character(object[[2]]),"/")[[3]],
                          dose = as.character(object[[3]]))
  data.e[,"row.num"]<-1:dim(data.e)[1]
  tmp.data <- list()
  for(i in 1:R){
    sampled.expand <- data.e[as.numeric(unlist(aggregate(row.num ~ data.e[,as.character(object[[3]])], data=data.e, 
                                                        FUN=function(x) sample(x,replace=TRUE))[[2]])),]
    tmp.data[[i]] <- aggregateBinomial(object, sampled.expand)
  }
    }
    if(type %in% c("continuous","Poisson","negbin1","negbin2")){
      data.e<-data
      data.e[,"row.num"]<-1:dim(data.e)[1]
      data.e[,"dose"]<-data.e[,as.character(object[[3]])]
      tmp.data <- list()
      for(i in 1:R){
        tmp.data[[i]] <- data.e[as.numeric(unlist(aggregate(row.num ~ dose, data=data.e, 
                                                             FUN=function(x) sample(x,replace=TRUE))[[2]])),]
         }
      }
  } else if(boot=="pseudorandom"){
    if(type=="binomial"){
    Y <- data[,strsplit(as.character(object[[2]]),"/")[[2]]]
    N <- data[,strsplit(as.character(object[[2]]),"/")[[3]]]
    shrinks <- which(Y==N | Y==0)
    Y[shrinks] <- Y[shrinks]+0.25
    N[shrinks] <- N[shrinks]+0.5
    prob <- rep(Y/N,N)
    tmp.data <- list()
    for(i in 1:R){
      sampled.expand <- data.frame(number = rbinom(length(prob),1,prob), 
                                   dose = rep(data[,as.character(object[[3]])],N), 
                                   total = 1)
      df <- aggregateBinomial(number/total~dose, sampled.expand)
      colnames(df) <- c(as.character(object[[3]]),
                                        strsplit(as.character(object[[2]]),"/")[[2]],
                                        strsplit(as.character(object[[2]]),"/")[[3]])
      tmp.data[[i]] <- df  
    }
    }
    if(type=="continuous"){
      mean.Y <- aggregate(data[,as.character(object[[2]])]~data[,as.character(object[[3]])],FUN=function(x) mean(x,na.rm=TRUE))[,2]
      sd.Y <- aggregate(data[,as.character(object[[2]])]~data[,as.character(object[[3]])],FUN=function(x) sd(x,na.rm=TRUE))[,2]
      Dose<- aggregate(data[,as.character(object[[2]])]~data[,as.character(object[[3]])],FUN=function(x) length(!is.na(x)))[,1]
      N.dose<- aggregate(data[,as.character(object[[2]])]~data[,as.character(object[[3]])],FUN=function(x) length(!is.na(x)))[,2]
      tmp.data <- list()
      for(i in 1:R){
        sampled <- data.frame(y = rnorm(sum(N.dose),mean=rep(mean.Y,N.dose),sd=rep(sd.Y,N.dose)), 
                                     dose = rep(Dose,N.dose))
        colnames(sampled) <- c(as.character(object[[2]]), as.character(object[[3]]))
        tmp.data[[i]] <- sampled
        }
    }
  }
    bmd.list <- lapply(tmp.data, function(x){
    bmdIso(object, data=x, type=type, bmr = bmr, backgType = backgType, backg=backg,def=def)})
  
    resMat <- matrix(NA,1,2)
    resMat[1,1] <- bmdIso(object, data=data, type=type, bmr = bmr, backgType = backgType, backg=backg,def=def) # quantile(unlist(bmd.list),0.5)
    resMat[1,2] <- quantile(unlist(bmd.list),0.05)
    colnames(resMat) <- c("BMD", "BMDL")
    rownames(resMat) <- c("")
    cat("\n\n")
    resMat    
}

