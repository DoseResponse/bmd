mjust<-function (modelList,expressions, dataused, level=0.95, seType="san") {
  require(sandwich, quietly = TRUE)
  require(car, quietly = TRUE)
  require(RLRsim, quietly = TRUE)
  require(Matrix, quietly = TRUE)
  deltab<-function (object, g, func = g, ...)
  {
    if (!is.character(g))
      stop("The argument 'g' must be a character string")
    if(inherits(object,c("lmerMod","lme"))){
      para<-fixef(object)
      coefVec<-fixef(object)
      para.names<-sapply(strsplit(names(coefVec), ":"), "[[", 1)
      para.names[1] <- gsub("\\(Intercept\\)", "Intercept", para.names[1])
      names(para)<-para.names
    }else{
      para <- coef(object)
      if(inherits(object,"lm")){
        para.names <- names(coef(object))
        para.names[1] <- gsub("\\(Intercept\\)", "Intercept", para.names[1])
        names(para) <- para.names
      }else{
        if(inherits(object,"drc")){
          coefVec<-coef(object)
          #para.names<-sapply(strsplit(names(coefVec), ":"), "[[", 1)
          para.names.tmp<-gsub(":", "_",names(coefVec))
          para.names<-sapply(strsplit(para.names.tmp, "_\\(Intercept\\)"), "[[", 1)
          names(para)<-para.names
        }else{
          para.names <- names(para)
        }}}
    g <- parse(text = g)
    q <- length(para)
    for (i in 1:q) {
      assign(names(para)[i], para[i])
    }
    gd <- rep(0,q)
    for (i in 1:q) {
      gd[i] <- eval(D(g, names(para)[i]))
    }
    gd
  }
  makeIIDdecomp <- function(modelObject,g)
  {
    if(inherits(modelObject, "lmerMod")){
      numObsUsed <- length(predict(modelObject))
      #data.tmp<-dataused
      #data.tmp[is.na(data.tmp)]<-10
      allRepUnits<-unique(dataused[,names(getME(modelObject,"flist"))])
      repUnitsUsed<-unique(unlist(getME(modelObject,"flist")))
      naRepUnits<-as.numeric(setdiff(allRepUnits, repUnitsUsed))
      numInd<-length(repUnitsUsed)
      beta<-getME(modelObject,"beta")
      X<-getME(modelObject,"X")
      Y<-getME(modelObject,"y")
      Z<-getME(modelObject,"Z")
      A<-getME(modelObject,"A")
      Sigma<-sigma(modelObject)
      R<-diag(Sigma^2,numObsUsed)
      V<-sigma(modelObject)^2*t(A)%*%A+R
      Vminus<-solve(V)
      f<-function(i){
        if(allRepUnits[i]%in%repUnitsUsed){
          tmpList<-which(getME(modelObject,"flist")[[1]]==allRepUnits[i])
          Xi<-matrix(X[tmpList,],nrow=length(tmpList))
          Zi<-matrix(Z[tmpList,],nrow=length(tmpList))
          Vi<-V[tmpList,tmpList]
          Yi<-matrix(Y[tmpList],nrow=length(tmpList))
          as.matrix(-t(Xi)%*%solve(Vi)%*%(Yi-Xi%*%beta))} else{matrix(rep(0,ncol(X)),ncol=1)}}
      EstFun<-matrix(unlist(lapply(seq_along(allRepUnits),f)),nrow=length(allRepUnits),byrow=T)
      db<-deltab(modelObject,g)
      iidVec0<-as.matrix(-db%*%solve(t(X)%*%Vminus%*%X)*length(repUnitsUsed))%*%t(EstFun)
      if (!is.null(naRepUnits)) {
        iidVec <- sqrt(length(allRepUnits)/numInd) * iidVec0
      } else {
        iidVec <- iidVec0
      }
    }else{
      if(inherits(modelObject, "lme")){
        numObsUsed <- length(predict(modelObject))
        allRepUnits<-unique(modelObject$data[,attr(getGroups(modelObject),"label")])
        repUnitsUsed<-unique(getGroups(modelObject))
        naRepUnits<-as.numeric(setdiff(allRepUnits, repUnitsUsed))
        numInd<-length(repUnitsUsed)
        beta<-fixed.effects(modelObject)
        X<-extract.lmeDesign(modelObject)$X
        Y<-as.matrix(extract.lmeDesign(modelObject)$y)
        Z<-extract.lmeDesign(modelObject)$Z
        GB<-getVarCov(modelObject)
        G.Block<-matrix(as.numeric(GB),nrow=sqrt(length(GB)))
        G<-bdiag(rep(list(G.Block),numInd))
        Sigma<-modelObject$sigma
        if(is.null(modelObject$modelStruct$corStruct)){
          R<-diag(Sigma^2,numObsUsed)
        }else{
          R<-Sigma^2*bdiag(corMatrix(modelObject$modelStruct$corStruct))
        }
        V<-Z%*%G%*%t(Z)+R
        Vminus<-solve(V)
        f<-function(i){
          if(allRepUnits[i]%in%repUnitsUsed){
            tmpList<-which(getGroups(modelObject)==allRepUnits[i])
            Xi<-matrix(X[tmpList,],nrow=length(tmpList))
            Zi<-matrix(Z[tmpList,],nrow=length(tmpList))
            Vi<-V[tmpList,tmpList]
            Yi<-matrix(Y[tmpList],nrow=length(tmpList))
            as.matrix(-t(Xi)%*%solve(Vi)%*%(Yi-Xi%*%beta))} else{matrix(rep(0,ncol(X)),ncol=1)}}
        EstFun<-matrix(unlist(lapply(seq_along(allRepUnits),f)),nrow=length(allRepUnits),byrow=T)
        db<-deltab(modelObject,g)
        iidVec0<-as.matrix(-db%*%solve(t(X)%*%Vminus%*%X)*length(repUnitsUsed))%*%t(EstFun)
        if (!is.null(naRepUnits)) {
          iidVec <- sqrt(length(allRepUnits)/numInd) * iidVec0
        } else {
          iidVec <- iidVec0
        }
      }else{
        numObsUsed <- ifelse(inherits(modelObject, "coxph"),
                             modelObject$n, ifelse(inherits(modelObject,"nls"),
                                                   length(predict(modelObject)),ifelse(inherits(modelObject,"drc"),
                                                                                       length(predict(modelObject)), nrow(modelObject$model))))
        db<-deltab(modelObject,g)
        iidVec0 <- db %*%bread(modelObject) %*% t(estfun(modelObject))
        moNAac <- modelObject$na.action
        numObs <- numObsUsed + length(moNAac)
        numInd<-numObs
        iidVec <- rep(0, numObs)
        if (!is.null(moNAac)) {
          iidVec[-moNAac] <- sqrt(numObs/numObsUsed) * iidVec0[!is.na(iidVec0)]
        } else {
          iidVec <- iidVec0
        }}}
    list(iidVec = iidVec, numObsUsed = numObsUsed, numInd = numInd)
  }
  iidList <- list()
  numModels <- length(modelList)
  for(i in 1:numModels)
  {
    iidList[[i]]<- makeIIDdecomp(modelList[[i]], expressions[[i]])
  }
  iidresp <- matrix(as.vector(unlist(lapply(iidList, function(listElt){listElt[[1]]}))), nrow = numModels, byrow = TRUE)
  numObsUsed <- as.vector(unlist(lapply(iidList, function(listElt) {
    listElt[[2]]})))
  thetaEst <- rep(NA, numModels)
  thetaSe <- rep(NA, numModels)
  for(i in 1:numModels)
  {
    if (inherits(modelList[[i]], "drc")){
      coefVec <- coef(modelList[[i]])
      #names(coefVec) <- sapply(strsplit(names(coefVec), ":"), "[[", 1)
      para.names.tmp<-gsub(":", "_",names(coefVec))
      names(coefVec)<-sapply(strsplit(para.names.tmp, "_\\(Intercept\\)"), "[[", 1)
      deltaRes <- deltaMethod(coefVec,expressions[[i]],vcov(modelList[[i]]))
    } else {
      if (inherits(modelList[[i]], "lmerMod")){
        coefVec <- fixef(modelList[[i]])
        names(coefVec) <- sapply(strsplit(names(coefVec), ":"), "[[", 1)
        deltaRes <- deltaMethod(coefVec,expressions[[i]],vcov(modelList[[i]]))
      }else{
        deltaRes <- deltaMethod(modelList[[i]],expressions[[i]])
      }}
    thetaEst[i] <- deltaRes[1]
    thetaSe[i] <- deltaRes[2]
  }
  thetaEst <- unlist(thetaEst)
  thetaSe <- unlist(thetaSe)
  ## Calculating the estimated variance-covariance matrix of the parameter estimates
  numInd <- iidList[[1]]$numInd
  covar <- (iidresp %*% t(iidresp)) / numInd
  vcMat <- covar / numInd # Defining the finite-sample variance-covariance matrix
  ## Replacing sandwich estimates by model-based standard errors
  modbas <- seType == "mod"
  if (any(modbas))
  {
    corMat <- cov2cor(vcMat)
    ## Retrieving standard errors for the specified estimate from the individual fits
    modSE <- thetaSe
    sanSE <- sqrt(diag(vcMat))
    sanSE[modbas] <- modSE[modbas]
    vcMat <- diag(sanSE,nrow=length(sanSE)) %*% corMat %*% diag(sanSE,nrow=length(sanSE))
  }
  quant <- qnorm(1 - (1 - level)/2)
  numInd <- iidList[[1]]$numInd
  varMA <- vcMat
  seMA <- sqrt(diag(varMA))
  quantVal <- quant * seMA
  zVec <- thetaEst*(1/seMA)
  pvals <- 1 - pchisq(zVec * zVec, 1)
  retMat <- as.matrix(cbind(thetaEst, seMA, thetaEst - quantVal, thetaEst + quantVal,pvals))
  colnames(retMat) <- c("Estimate", "Std. Error", "Lower", "Upper", "Pr(>|z|)")
  output<-list(retMat,varMA)
  names(output)<-c("coef","covar")
  return(invisible(output))
}
