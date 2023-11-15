bmdOrdinal <- function(object, bmr=0.1, backgType = "modelBased", def="excess", interval = "delta", level = 0.95, R = 500, bootType = "nonparametric", display = TRUE){
  bmdList <- lapply(object$drmList, function(mod) bmd(mod, bmr = bmr, backgType = backgType, def=def, display=FALSE))
  BMD <- mean(sapply(bmdList, FUN=function(x) x$Results[1]))
  
  if(interval == "delta"){  
    if(substr(object$drmList[[1]]$fct$name, 1,2) %in% c("FP")){
      #cat("parametric CI not available for LN and FPL models.\n")
      warning("parametric CI not available for LN and FPL models.")
      CI <- NA
    } else {
      if(object$drmList[[1]]$fct$name == "LL.2"){
        BMD.pooled <- mjust(object$drmList,
                            as.list(rep(paste0("exp(log(1/",bmr,"-1)/b+log(e))"), length(object$drmList))), # closed solution for LL.2
                            seType = "san")
      } else if(object$drmList[[1]]$fct$name == "LN.2"){
        constant <- qnorm(bmr)
        BMD.pooled <- mjust(object$drmList,
                            as.list(rep(paste0("exp(",constant,"/b)*e"), length(object$drmList))), # closed solution for LN.2
                            seType = "san")
      } else if(object$drmList[[1]]$fct$name == "W1.2"){
        BMD.pooled <- mjust(object$drmList,
                            as.list(rep(paste0("(-log(",bmr,"))^(1/b)*e"), length(object$drmList))), # closed solution for W1.2
                            seType = "san")
      } else if(object$drmList[[1]]$fct$name == "W2.2"){
        BMD.pooled <- mjust(object$drmList,
                            as.list(rep(paste0("(-log(1-",bmr,"))^(1/b)*e"), length(object$drmList))), # closed solution for W2.2
                            seType = "san")
      }
      tmp <- confint(multcomp:::glht(multcomp:::parm(BMD.pooled[["coef"]][,1],
                                                     BMD.pooled[["covar"]]), linfct = matrix(rep(1/length(object$drmList), length(object$drmList)),1,length(object$drmList))),
                     level = level)
      CI <- tmp$confint[1,2:3]
    }
  } else if(interval == "bootstrap") {
    bootData <- bootDataGenOrdinal(object, R=R, bootType = bootType)
    
    bmdBoot <- numeric(R)
    for(i in 1:R){
      modelBoot <- suppressWarnings(try(drmOrdinal(object$levels, object$dose, object$weights, bootData[[i]], object$fct), silent = TRUE))
      bmdAllBoot <- lapply(modelBoot$drmList, function(mod) try(bmd(mod, bmr = bmr, backgType = backgType, def=def, display=FALSE)$Results[1], silent = TRUE))
      bmdBoot[i] <- mean(as.numeric(bmdAllBoot))
    }
    CI <- quantile(bmdBoot, c(1-level, 1 - level), na.rm = TRUE)
  } else {
    CI <- NA
  }
  
  resMat <- matrix(data = c(BMD, CI[1]), ncol = 2, dimnames=list(NULL, c("BMD", "BMDL")))
  
  resBMD<-list(Results = resMat,
               interval = CI,
               bmdList = bmdList)
  
  class(resBMD) <- "bmdOrdinal"
  if(display){ print(resMat) }
  invisible(resBMD)
}
