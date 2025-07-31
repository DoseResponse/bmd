invBmd <- function(object, bmr, level=0.9, slope, backgType="modelBased", 
                     backg=NA, catLev=NA, extFactor=10, def, useSD=useSD, sandwich.vcov=FALSE){
  
  ParmVec0 <- object$fct$fixed
  ParmVec <- ParmVec0
  ParmVec[is.na(ParmVec0)] <- coef(object)
  g<-object$fct$fct
  h<-object$fct$deriv1
  if(substr(deparse(as.list(body(h))[[length(as.list(body(h)))]])[1],1,6)=="return"){
    element.num <- length(as.list(body(h)))-1
    dpgh <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
    pgh <- paste(dpgh, "[, notFixed]", sep="")
    body(h)[[element.num]] <- as.call(str2lang(pgh))
    body(h)[[element.num+1]] <- substitute(derMat)
  } else {
    element.num <- length(as.list(body(h)))
  }
  
  
  if(def %in% c("additional", "added")){
    if(identical(slope,"decreasing")){
      g0 <- function(x){as.numeric(g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                     g(x, matrix(coef(object), 1, length(coef(object)))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " + as.matrix(cbind(0, 0, -1, 0))[, notFixed]", sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " + as.matrix(cbind(0, 0, -1, 0, 0))[, notFixed]", sep="")
      }
      body(h)[[element.num]] <- as.call(str2lang(pgh))
      
    } else if(identical(slope,"increasing")){
      g0 <- function(x){as.numeric(g(x, matrix(coef(object), 1, length(coef(object)))) - 
                                     g(0, matrix(coef(object), 1, length(coef(object)))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " + as.matrix(cbind(0, -1, 0, 0))[, notFixed]", sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " + as.matrix(cbind(0, -1, 0, 0, 0))[, notFixed]", sep="")
      }
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  if(def == "excess"){
    if(identical(slope,"decreasing")){
      g0 <- function(x){as.numeric((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(x, matrix(coef(object), 1, length(coef(object)))))/
                                     g(0, matrix(coef(object), 1, length(coef(object)))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(-1/(parmMat[,3]), 
                   -1/(1-parmMat[,3]), 
                   -1/(1-parmMat[,3]), -1/(1-parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 0, (", 
                     dpg, ") /(parmMat[,3]^2),
                   0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(-1/(parmMat[,3]), 
                   -1/(1-parmMat[,3]), -1/(1-parmMat[,3]), 
                   -1/(1-parmMat[,3]), -1/(1-parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 0, (", 
                     dpg, ") /(parmMat[,3]^2),
                   0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    } else if(identical(slope,"increasing")){
      g0 <- function(x){as.numeric((g(x, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(0, matrix(coef(object), 1, length(coef(object)))))/
                                     (1 - g(0, matrix(coef(object), 1, length(coef(object))))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(1-parmMat[,2]), 
                   1/(1-parmMat[,2]), 
                   1/(1-parmMat[,2]), 1/(1-parmMat[,2])))[, notFixed]",
                     " + as.matrix(cbind(0, (", 
                     dpg, " - 1)/(1-parmMat[,2])^2,
                   0, 0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(1-parmMat[,2]), 
                   1/(1-parmMat[,2]), 
                   1/(1-parmMat[,2]), 1/(1-parmMat[,2]), 1/(1-parmMat[,2])))[, notFixed]",
                     " + as.matrix(cbind(0, (", 
                     dpg, " - 1)/(1-parmMat[,2])^2,
                   0, 0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  if(def == "relative"){
    if(identical(slope,"decreasing")){
      g0 <- function(x){as.numeric((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(x, matrix(coef(object), 1, length(coef(object)))))/
                                     g(0, matrix(coef(object), 1, length(coef(object)))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]", dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(-1/(parmMat[,3]), 
                   -1/(parmMat[,3]), 
                   -1/(parmMat[,3]), -1/(parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 0, (", 
                     dpg, ") /(parmMat[,3]^2),
                   0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(-1/(parmMat[,3]), 
                   -1/(parmMat[,3]), -1/(parmMat[,3]), 
                   -1/(parmMat[,3]), -1/(parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 0, (", 
                     dpg, ") /(parmMat[,3]^2),
                   0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
      
    } else if(identical(slope,"increasing")){
      g0 <- function(x){as.numeric((g(x, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(0, matrix(coef(object), 1, length(coef(object)))))/
                                     (g(0, matrix(coef(object), 1, length(coef(object))))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,2]), 
                   1/(parmMat[,2]), 
                   1/(parmMat[,2]), 1/(parmMat[,2])))[, notFixed]",
                     " - as.matrix(cbind(0, (", 
                     dpg, " )/(parmMat[,2]^2),
                   0, 0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,2]), 
                   1/(parmMat[,2]), 
                   1/(parmMat[,2]), 1/(parmMat[,2]), 1/(parmMat[,2])))[, notFixed]",
                     " - as.matrix(cbind(0, (", 
                     dpg, ")/(parmMat[,2]^2),
                   0, 0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  if(def == "point"){
    stop(paste("Inverse regression not possible for def=point", sep=""))
  }
  
  if(def == "extra"){
    if(identical(slope,"decreasing")){
      g0 <- function(x){as.numeric((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(x, matrix(coef(object), 1, length(coef(object)))))/
                                     (g(0, matrix(coef(object), 1, length(coef(object)))) 
                                      - g(Inf, matrix(coef(object), 1, length(coef(object))))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,2]-parmMat[,3]), 
                   1/(parmMat[,2]-parmMat[,3]), 
                   1/(parmMat[,2]-parmMat[,3]), 1/(parmMat[,2]-parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (", dpg, " - parmMat[,2])/(parmMat[,2]-parmMat[,3])^2,
                     ( parmMat[,3] - (", dpg, "))/(parmMat[,2]-parmMat[,3])^2,
                     0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,2]-parmMat[,3]), 
                   1/(parmMat[,2]-parmMat[,3]), 1/(parmMat[,2]-parmMat[,3]), 
                   1/(parmMat[,2]-parmMat[,3]), 1/(parmMat[,2]-parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (", dpg, " - parmMat[,2])/(parmMat[,2]-parmMat[,3])^2,
                     ( parmMat[,3] - (", dpg, "))/(parmMat[,2]-parmMat[,3])^2,
                     0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    } else if(identical(slope,"increasing")){
      g0 <- function(x){as.numeric((g(x, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(0, matrix(coef(object), 1, length(coef(object)))))/
                                     (g(Inf, matrix(coef(object), 1, length(coef(object)))) 
                                      - g(0, matrix(coef(object), 1, length(coef(object))))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,3]-parmMat[,2]), 
                   1/(parmMat[,3]-parmMat[,2]), 
                   1/(parmMat[,3]-parmMat[,2]), 1/(parmMat[,3]-parmMat[,2])))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (", dpg, " - parmMat[,3])/(parmMat[,3]-parmMat[,2])^2,
                     ( parmMat[,2] - (", dpg, "))/(parmMat[,3]-parmMat[,2])^2,
                     0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,3]-parmMat[,2]), 
                   1/(parmMat[,3]-parmMat[,2]), 
                   1/(parmMat[,3]-parmMat[,2]), 1/(parmMat[,3]-parmMat[,2]), 1/(parmMat[,3]-parmMat[,2])))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (", dpg, " - parmMat[,3])/(parmMat[,3]-parmMat[,2])^2,
                     ( parmMat[,2] - (", dpg, "))/(parmMat[,3]-parmMat[,2])^2,
                     0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  ###
  if((def == "hybridAdd" | def == "hybridExc") & !(backgType %in% c("hybridSD","absolute"))){
    stop(paste("When def = ", def, "backgType needs to be either hybridSD or absolute"))
  }
  
  if(def == "hybridAdd" & identical(backgType,"hybridSD")){
    if(identical(slope,"decreasing")){
      NumSD <- ifelse(is.na(backg), 2, backg)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric(pnorm((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                            NumSD*useSD-g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - 
                                     p0)}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," - (",dpg,"))/(",useSD,")))/",useSD,", 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,", 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    } else if(identical(slope,"increasing")){
      NumSD <- ifelse(is.na(backg), 2, backg)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric(1 - pnorm((NumSD*useSD + g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                                g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - p0)}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]", dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,", 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,", 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  if(def == "hybridAdd" & identical(backgType,"absolute")){
    if(is.na(backg)){
      stop(paste("No value of backg is provided", sep=""))
    }
    if(identical(slope,"decreasing")){
      NumSD <- abs((backg - g(0, matrix(coef(object), 1, length(coef(object)))))/useSD)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric(pnorm((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                            NumSD*useSD-g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - 
                                     p0)}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," - (",dpg,"))/(",useSD,")))/",useSD,", 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,", 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    } else if(identical(slope,"increasing")){
      NumSD <- abs((backg - g(0, matrix(coef(object), 1, length(coef(object)))))/useSD)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric(1 - pnorm((NumSD*useSD + g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                                g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - p0)}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]", dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,", 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,", 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  
  ###
  
  if(def == "hybridExc" & identical(backgType,"hybridSD")){
    if(identical(slope,"decreasing")){
      NumSD <- ifelse(is.na(backg), 2, backg)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric((pnorm((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                             NumSD*useSD-g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - 
                                      p0)/(1-p0))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," - (",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/(",useSD,"*(1-",p0,")))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/(",useSD,"*(1-",p0,")))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
      
    } else if(identical(slope,"increasing")){
      NumSD <- ifelse(is.na(backg), 2, backg)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric((1 - pnorm((NumSD*useSD + g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                                 g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - p0)/(1-p0))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]", dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,"))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,"))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  if(def == "hybridExc" & identical(backgType,"absolute")){
    if(is.na(backg)){
      stop(paste("No value of backg is provided", sep=""))
    }
    if(identical(slope,"decreasing")){
      NumSD <- abs((backg - g(0, matrix(coef(object), 1, length(coef(object)))))/useSD)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric((pnorm((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                             NumSD*useSD-g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - 
                                      p0)/(1-p0))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," - (",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/(",useSD,"*(1-",p0,")))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/(",useSD,"*(1-",p0,")))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
      
    } else if(identical(slope,"increasing")){
      NumSD <- abs((backg - g(0, matrix(coef(object), 1, length(coef(object)))))/useSD)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric((1 - pnorm((NumSD*useSD + g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                                 g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - p0))/(1-p0)}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]", dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,"))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,"))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  ###  
  respType <- object[["type"]]
  if (identical(respType, "continuous"))
  {
    tquan <- qt(1 - (1 - level)/2, df.residual(object))   
  } else {
    tquan <- qnorm(1 - (1 - level)/2)
  }
  
  
  
  #  g1<-function(x){as.numeric(g(x, matrix(coef(object), 1, length(coef(object)))))}
  if(sandwich.vcov){
    j1<-function(x){sqrt(as.vector(h(x, matrix(coef(object), 1, length(coef(object)))))%*%
                           sandwich::sandwich(object)%*%
                           as.vector(h(x, matrix(coef(object), 1, length(coef(object))))))}
  } else{
    j1<-function(x){sqrt(as.vector(h(x, matrix(coef(object), 1, length(coef(object)))))%*%
                           vcov(object)%*%
                           as.vector(h(x, matrix(coef(object), 1, length(coef(object))))))}
  }
  my.fun1<-function(x){g0(x) + tquan * j1(x)}
  my.fun2<-function(x){g0(x) - tquan * j1(x)}
  
  objDL <- object[["dataList"]][["names"]]
  #colnames(newData0) <- c(objDL[["dName"]], objDL[["cName"]])
  maxdose <- extFactor * max(object[["dataList"]][["dose"]])
  
  
  rootFctBMD <- function(x) 
  {
    newData <- data.frame(x, catLev)
    colnames(newData) <- c(objDL[["dName"]], objDL[["cName"]])
    #    print(c(x, predict(object, newData, interval = intType, level = level)[2] - yval))
    g0(x) - bmr
  }
  
  rootFct1 <- function(x) 
  {
    newData <- data.frame(x, catLev)
    colnames(newData) <- c(objDL[["dName"]], objDL[["cName"]])
    my.fun1(x) - bmr
  }
  
  rootFct2 <- function(x) 
  {
    newData <- data.frame(x, catLev)
    colnames(newData) <- c(objDL[["dName"]], objDL[["cName"]])
    my.fun2(x) - bmr
  }
  
  ResMat <- matrix(NA,1,3)
  colnames(ResMat)<-c("BMD","BMDL","BMDU")
  BMDLevel<-try(uniroot(rootFctBMD, c(0, maxdose)), silent = TRUE)
  if(inherits(BMDLevel,"try-error")){
    stop(paste("No solution possible for this BMR", sep=""))
  } else {
    ResMat[1,1] <- BMDLevel$root
    
    upLevel<-try(uniroot(rootFct1, c(0, maxdose)), silent = TRUE)
    ResMat[1,2]<- ifelse(inherits(upLevel,"try-error"), 0, upLevel$root)
    
    lowLevel<-try(uniroot(rootFct2, c(0, maxdose)), silent = TRUE)
    ResMat[1,3]<- ifelse(inherits(lowLevel,"try-error"), Inf, lowLevel$root)
  }
  return(ResMat)
}

