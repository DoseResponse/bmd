bmdIso <- function(object, data, type, bmr, p0, backgType = c("modelBased", "absolute","hybridSD","hybridPercentile"), backg=NA, def = c("excess", "additional", "relative", "added", "hybridExc", "hybridAdd", "point"), display=FALSE){
  PAV.p <- PAV(object, data, type)
  n <- as.numeric(table(data[, paste(object[[3]])]))
  if(type=="continuous"){
    #sigma.sq <- sum((data[,paste(object[[2]])]-rep(PAV.p,n))^2)/length(data[,paste(object[[2]])])
    sigma.sq <- sd(data[data[,paste(object[[3]])]==0,paste(object[[2]])])
  }
  Dose <- sort(unique(data[,paste(object[[3]])]))
  f0 <- min(PAV.p) 
  if (missing(def)) {
    stop(paste("def is missing", sep=""))
  }
  if (!(def %in% c("excess", "additional", "relative", "added", "hybridExc", "hybridAdd", "point"))) {
    stop(paste("Could not recognize def", sep=""))
  }
  if (identical(backgType,"modelBased")) {
    background <- f0
  } else if (identical(backgType,"absolute") & !(def %in% c("relative"))) {
    background <- backg
  } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
    background <- ifelse(is.na(backg),0,backg)
  } else if (identical(backgType,"hybridSD")) {
    background <- ifelse(is.na(backg), 1-pnorm(2), 1-pnorm(backg))
  } else if (identical(backgType,"absolute") & 
             (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
    background <- ifelse(is.na(backg), 
                         1-pnorm(2),
                         1-pnorm((backg-f0)/sigma.sq))
  } else {
    background <- ifelse(is.na(backg),1-0.9,1-backg)
  }
  def <- match.arg(def)
  if (identical(type, "binomial")) {
    Cq <- switch(def, 
                 excess = bmr * (1 - background) + background, 
                 additional = bmr + background, 
                 point = bmr)
  }
  if (identical(type, "binomial") & (def %in% c("relative","added", "hybridExc","hybridAdd"))) {
    stop(paste("\"",def, "\" is not available for quantal data", sep=""))
  }
  if (type %in% c("Poisson","negbin1","negbin2")) {
      Cq <- switch(def,
                   relative = bmr * background + background,
                   extra = bmr*abs(diff(predict(object, data.frame(c(0, Inf))))) + background,
                   point = bmr)
    }
  if (type %in% c("Poisson","negbin1","negbin2") & (def %in% c("excess","additional","added","hybridExc","hybridAdd"))) {
    stop(paste("\"",def, "\" is not available for count data", sep=""))
  }
  if (identical(type, "continuous") & (def %in% c("excess", "additional"))) {
    stop(paste("\"",def, "\" is not available for continuous data", sep=""))
  }
  if (identical(type, "continuous")) {
    Cq <- switch(def, relative = bmr * background + background, 
                 added = bmr + background,
                 point = bmr,
                 hybridAdd = sigma.sq * 
                   (qnorm(1 - background) - qnorm(1 - (background + bmr))) + 
                   f0,
                 hybridExc = sigma.sq * 
                   (qnorm(1 - background) - qnorm(1 + background - (1 - background)*bmr)) + 
                   f0)
  } 
  if(sum(Cq == PAV.p)>0){
    BMD <- max(Dose[Cq == PAV.p])
  } else if( length(PAV.p) >  sum(as.numeric(Cq > PAV.p)) ){
    crit <- sum(as.numeric(Cq > PAV.p))
    BMD <- Dose[crit] + (Cq-PAV.p[crit])*(Dose[crit+1]-Dose[crit])/(PAV.p[crit+1]-PAV.p[crit])  
  } else if(sum(as.numeric(Cq < PAV.p))==0){
    BMD <- 0
  } else if(sum(as.numeric(Cq < PAV.p)) == length(PAV.p)){
    BMD <- Dose[length(PAV.p)]
  }
  if (display) {
    cat("Effective response level: ", Cq )
  }
  BMD
}
