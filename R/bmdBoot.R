bmdBoot <- function(object, bmr, R=1000, boot="nonparametric", bmdType = "orig",
                    backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
                    backg=NA, 
                    def = c("excess", "additional", 
                            "relative", "extra", "added", "hybridExc", "hybridAdd", "point")){
  tmp.data <- bootDataGen(object,R,boot)
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


