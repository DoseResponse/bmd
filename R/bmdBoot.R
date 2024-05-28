bmdBoot <- function(object, bmr, R=1000, bootType="nonparametric", bmdType = "orig",
                    backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
                    backg=NA, 
                    controlSD=NA,
                    def = c("excess", "additional", 
                            "relative", "extra", "added", "hybridExc", "hybridAdd", "point"),
                    respTrans = c("none", "log", "sqrt"),
                    bootInterval = c("percentile","BCa"),
                    display=TRUE, level=0.95){
  if (identical(object$type,"binomial") & bootType=="semiparametric") {
    stop(paste("\"Semiparametric bootstrap does not work for quantal data\"", sep=""))
  }
  if (object$type %in% c("Poisson","negbin1","negbin2") & bootType!="nonparametric") {
    stop(paste("\"",object$type,"\" only works with nonparametric bootstrap\"", sep=""))
  }
  
  get.drm.list <- function(tmp.data){
    if(ncol(object$parmMat) == 1){
      drm.list.tmp <- lapply(tmp.data, function(x){
        try(drm(object$call$formula, data = x, type = object$type, fct = object[["fct"]]), TRUE)
      }
      )
    } else if(is.null(object$call$pmodels)){
      drm.list.tmp <- lapply(tmp.data, function(x){
        if(object$type != "binomial"){
          x[[as.character(object$call$curveid)]] <- x[[paste0("orig.", as.character(object$call$curveid))]]
        }
        try(
          eval(substitute(drm(object$call$formula, weights = weights0, curveid = curveid0,
                              data = x, type = object$type, fct = object$fct, control = drmc(noMessage = TRUE)),
                          list(weights0 = object$call$weights,
                               curveid0 = object$call$curveid)
          )),
          silent = TRUE)
      })
    } else {
      drm.list.tmp <- lapply(tmp.data, function(x){
        if(object$type != "binomial"){
          x[[as.character(object$call$curveid)]] <- x[[paste0("orig.", as.character(object$call$curveid))]]
        }
        try(
          eval(substitute(drm(object$call$formula, weights = weights0, curveid = curveid0,pmodels = pmodels0,
                              data = x, type = object$type, fct = object$fct, control = drmc(noMessage = TRUE)),
                          list(weights0 = object$call$weights,
                               curveid0 = object$call$curveid,
                               pmodels0 = object$call$pmodels)
          )),
          silent = TRUE)
      })
    }
  }
  
  get.bmd <- function(x){
    if(ncol(object$parmMat) == 1){
      try(bmd(x, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
              display=FALSE, level=level)[["Results"]][1],TRUE)
    } else {
      try(bmd(x, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
              display=FALSE, level=level)[["Results"]][colnames(object$parmMat), 1],TRUE)
    }
  }
 
  if (object$type %in% c("binomial","continuous")) {
    
    tmp.data <- bootDataGen(object,R,bootType,aggregated=FALSE)
    
    drm.list.tmp <- get.drm.list(tmp.data)
    list.condition <- sapply(drm.list.tmp, function(x) class(x)=="drc")
    drm.list  <- drm.list.tmp[list.condition]
    
    bmd.list.try <- lapply(drm.list,get.bmd)
    bmd.list <- bmd.list.try[!sapply(bmd.list.try, function(x) any(is.na(x)))]
  }
  
  if (object$type %in% c("Poisson","negbin1","negbin2")) {
    tmp.data <- bootDataGen(object,R,bootType,aggregated=FALSE)
  
    # drm.list.tmp <- lapply(tmp.data, function(x){
    #   try(drm(object$call$formula, data = x, type = object$type, weights=weights, fct = object[["fct"]]),TRUE)}
    # )
    drm.list.tmp <- get.drm.list(tmp.data)
    list.condition <- sapply(drm.list.tmp, function(x) class(x)=="drc")
    drm.list  <- drm.list.tmp[list.condition]
    
    # bmd.list <- lapply(drm.list,function(x){
    #   bmd(x, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
    #       display=FALSE, level=level)[["Results"]][1]}
    # )
    bmd.list.try <- lapply(drm.list,get.bmd)
    bmd.list <- bmd.list.try[!sapply(bmd.list.try, function(x) any(is.na(x)))]
  }

  
  if(identical(object$type, "continuous")){
    if(identical(bootInterval,"BCa")){
      jackData <- list()
      for(i in 1:(dim(object$data)[1])){
        jackData[[i]] <- object$data[-i,]
      }
      # bootJack.drm.tmp <- lapply(jackData, function(x){
      #   try(drm(object$call$formula, data = x, fct = object[["fct"]]),TRUE)
      # })
      bootJack.drm.tmp <- get.drm.list(jackData)
      list.condition <- sapply(bootJack.drm.tmp, function(x) class(x)=="drc")
      bootJack.drm<- bootJack.drm.tmp[list.condition]
      
      # bootJack <- sapply(bootJack.drm, function(x){
      #   bmd(x, bmr, backgType = backgType, backg = backg, def = def, controlSD=controlSD, respTrans = respTrans,
      #       interval = "delta", display=FALSE, level=level)$Results[1]
      # }
      # )
      
      bootJack.list.try <- lapply(bootJack.drm, get.bmd)
      bootJack.list <- bootJack.list.try[!sapply(bootJack.list.try, function(x) any(is.na(x)))]
      
      # use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
      #                display=FALSE, level=level)[["Results"]][1]
      use.bmd <- get.bmd(object)
      if(ncol(object$parmMat)==1){
        BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = object$data, unlist(bmd.list), unlist(bootJack.list), level=level)[1])
      } else {
        BCaBMDL <- sapply(1:ncol(object$parmMat), function(i) as.numeric(BCa(obs = use.bmd[i], data = object$data, 
                                                                             sapply(bmd.list, function(x) x[i]), 
                                                                             sapply(bootJack.list, function(x) x[i]), level = level)[1]))
      }
      
    }
  }
  if(identical(object$type, "binomial")){
    if(identical(bootInterval, "BCa")){
      # data.str <- object$data
      # data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
      # data.e<-expandBinomial(data.str, 
      #                        number = "number",
      #                        total = "weights",
      #                        dose = as.character(object$call$formula[[3]]))
      
      if(ncol(object$parmMat) == 1){
        data.str <- object$data
        data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
        data.e<-expandBinomial(data.str, 
                               number = "number",
                               total = "weights",
                               dose = as.character(object$call$formula[[3]]))
        df <- data.frame(data.e[,as.character(object$call$formula[[3]])],
                         data.e[,"number"],
                         data.e[,"weights"])
        colnames(df) <- c(as.character(object$call$formula[[3]]),
                          as.character(object$call$formula[[2]])[[2]],
                          as.character(object$call$formula[[2]])[[3]])
        df$rowNum <-1:nrow(df)
      } else {
        data.str <- object$data
        data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
        data.e<-expandBinomial(data.str, 
                               number = "number",
                               total = "weights",
                               dose = as.character(object$call$formula[[3]]),
                               curveid = as.character(object$call$curveid))
        df <- data.frame(data.e[,as.character(object$call$formula[[3]])],
                         data.e[,"number"],
                         data.e[,"weights"],
                         data.e[,as.character(object$call$curveid)])
        colnames(df) <- c(as.character(object$call$formula[[3]]),
                          as.character(object$call$formula[[2]])[[2]],
                          as.character(object$call$formula[[2]])[[3]],
                          as.character(object$call$curveid))
      }
      
      jackData <- list()
      for(i in 1:(dim(data.e)[1])){
        jackData[[i]] <- data.e[-i,]
      }
      
      bootJack.drm.tmp <- lapply(jackData, function(x){
        try(drm(number~dose, data = x, type = "binomial", fct = object[["fct"]]),TRUE)
      })
      bootJack.drm.tmp <- get.drm.list(jackData)
      list.condition <- sapply(bootJack.drm.tmp, function(x) class(x)=="drc")
      bootJack.drm<- bootJack.drm.tmp[list.condition]
      
      # bootJack <- sapply(bootJack.drm, function(x){
      #   bmd(x, bmr, backgType = backgType, backg = backg, def = def, controlSD=controlSD, respTrans = respTrans,
      #       interval = "delta", display=FALSE, level=level)$Results[1]
      # }
      # )
      
      bootJack.list.try <- lapply(bootJack.drm, get.bmd)
      bootJack.list <- bootJack.list.try[!sapply(bootJack.list.try, function(x) any(is.na(x)))]
      
      # use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
      #                display=FALSE, level=level)[["Results"]][1]
      # BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = data.e, unlist(bmd.list), bootJack, level=level)[1])
      
      use.bmd <- get.bmd(object)
      if(ncol(object$parmMat)==1){
        BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = object$data, unlist(bmd.list), unlist(bootJack.list), level=level)[1])
      } else {
        BCaBMDL <- sapply(1:ncol(object$parmMat), function(i) as.numeric(BCa(obs = use.bmd[i], data = object$data, 
                                                                             sapply(bmd.list, function(x) x[i]), 
                                                                             sapply(bootJack.list, function(x) x[i]), level = level)[1]))
      }
      
    }
  }
  if(object$type %in% c("Poisson","negbin1","negbin2")){
    if(identical(bootInterval,"BCa")){
      jackData <- list()
      for(i in 1:(dim(object$data)[1])){
        jackData[[i]] <- object$data[-i,]
      }
      # bootJack.drm.tmp <- lapply(jackData, function(x){
      #   try(drm(object$call$formula, data = x, type = object$type, weights=weights, fct = object[["fct"]]),TRUE)
      # })
      bootJack.drm.tmp <- get.drm.list(jackData)
      list.condition <- sapply(bootJack.drm.tmp, function(x) class(x)=="drc")
      bootJack.drm<- bootJack.drm.tmp[list.condition]
      
      # bootJack <- sapply(bootJack.drm, function(x){
      #   bmd(x, bmr, backgType = backgType, backg = backg, def = def, controlSD=controlSD, interval = "delta", respTrans = respTrans,
      #       display=FALSE, level=level)$Results[1]
      # }
      # )
      
      bootJack.list.try <- lapply(bootJack.drm, get.bmd)
      bootJack.list <- bootJack.list.try[!sapply(bootJack.list.try, function(x) any(is.na(x)))]
      
      # use.bmd <- bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans,
      #                display=FALSE, level=level)[["Results"]][1]
      # BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = object$data, bootSample=unlist(bmd.list), 
      #                           bootjack=bootJack, level=level)[1])
      use.bmd <- get.bmd(object)
      if(ncol(object$parmMat)==1){
        BCaBMDL <- as.numeric(BCa(obs = use.bmd, data = object$data, unlist(bmd.list), unlist(bootJack.list), level=level)[1])
      } else {
        BCaBMDL <- sapply(1:ncol(object$parmMat), function(i) as.numeric(BCa(obs = use.bmd[i], data = object$data, 
                                                                             sapply(bmd.list, function(x) x[i]), 
                                                                             sapply(bootJack.list, function(x) x[i]), level = level)[1]))
      }
    }
  }
  if(bmdType == "orig"){
    use.bmd <- get.bmd(object) # bmd(object, bmr = bmr, backgType = backgType, backg=backg, def=def, controlSD=controlSD, respTrans = respTrans, display=FALSE, level=level)[["Results"]][,1]
  } else if(bmdType == "mean"){
    if(ncol(object$parmMat) == 1){
      use.bmd <- mean(unlist(bmd.list))
    } else {
      use.bmd <- apply(matrix(unlist(bmd.list), ncol = ncol(object$parmMat), byrow = TRUE), 2, mean, na.rm = TRUE)
    }
  } else if(bmdType == "median"){
    if(ncol(object$parmMat) == 1){
      use.bmd <- quantile(unlist(bmd.list),0.5)  
    } else {
      use.bmd <- apply(matrix(unlist(bmd.list), ncol = ncol(object$parmMat), byrow = TRUE), 2, quantile, p= 0.5, na.rm = TRUE)
    }
  } 
  
  if(ncol(object$parmMat) == 1){
    BMDL <- ifelse(identical(bootInterval,"BCa"), BCaBMDL, quantile(unlist(bmd.list),c(1-level),na.rm = TRUE))
    BMDU <- ifelse(identical(bootInterval,"BCa"), "Not available for BCa bootstrap", quantile(unlist(bmd.list),c(level),na.rm = TRUE))
  } else {
    if(identical(bootInterval, "BCa")){
      BMDL <- BCaBMDL
      BMDU <- rep("Not available for BCa bootstrap", ncol(object$parmMat))
    } else {
      BMDL <- apply(matrix(unlist(bmd.list), ncol = ncol(object$parmMat), byrow = TRUE), 
                    2, quantile, p= c(1-level), na.rm = TRUE)
      BMDU <- apply(matrix(unlist(bmd.list), ncol = ncol(object$parmMat), byrow = TRUE), 2, quantile, p= c(level), na.rm = TRUE)
    }
  }
  
  resMat <- matrix(NA,ncol(object$parmMat),2)
  resMat[,1] <- use.bmd
  resMat[,2] <- BMDL
  colnames(resMat) <- c("BMD", "BMDL")
  
  intMat <- matrix(NA,ncol(object$parmMat),2)
  intMat[,1] <- BMDL
  intMat[,2] <- BMDU
  colnames(intMat) <- c("Lower", "Upper")
  
  if(ncol(object$parmMat) == 1){
    rownames(resMat) <- c("")
    rownames(intMat) <- c("")
  } else {
    rownames(resMat) <- colnames(object$parmMat)
    rownames(intMat) <- colnames(object$parmMat)
  }
  
  if(display){
    print(resMat)
  }
  
  resBMD<-list(Results = resMat,
               bootEst = unlist(bmd.list),
               Interval = intMat)
  class(resBMD) <- "bmd"
  invisible(resBMD)
}