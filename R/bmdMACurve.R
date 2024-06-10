bmdMACurve<-function(modelList,modelWeights,bmrScaled0, searchInterval="dataBased"){
  nCurves <- ncol(modelList[[1]]$parmMat)
  
  if(nCurves == 1){
    fList <- mtList()
    
    f.all<-list()
    fi<-lapply(modelList,function(x){fList[[x[["fct"]][["name"]]]]$fct})
    
    for(i in 1:length(modelList)){
      if(!identical(modelList[[i]]$fct$text,"Fractional polynomial")){
      if(length(modelList[[i]]$fct$fixed)==5){
        parm<-modelList[[i]]$fct$fixed
        parm[is.na(parm)]<-coef(modelList[[i]])
        body(fi[[i]])<-gsub("b1",parm[1],paste(body(fi[[i]])))[2]
        body(fi[[i]])<-gsub("c1",parm[2],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("d1",parm[3],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("e1",parm[4],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("f1",parm[5],paste(body(fi[[i]])))
        f.all[[i]] <- fi[[i]]
      } else if(length(modelList[[i]]$fct$fixed)==4){
        parm<-modelList[[i]]$fct$fixed
        parm[is.na(parm)]<-coef(modelList[[i]])
        body(fi[[i]])<-gsub("b1",parm[1],paste(body(fi[[i]])))[2]
        body(fi[[i]])<-gsub("c1",parm[2],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("d1",parm[3],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("e1",parm[4],paste(body(fi[[i]])))
        f.all[[i]] <- fi[[i]]
      } 
    } else if(identical(modelList[[i]]$fct$text,"Fractional polynomial")){
        fi[[i]] <- fList[["FPL.4"]]$fct
        parm<-modelList[[i]]$fct$fixed
        parm[is.na(parm)]<-coef(modelList[[i]])
        body(fi[[i]])<-gsub("b1",parm[1],paste(body(fi[[i]])))[2]
        body(fi[[i]])<-gsub("c1",parm[2],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("d1",parm[3],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("e1",parm[4],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("p1",unlist(strsplit(modelList[[i]]$fct$name, "\\,|\\(|\\)"))[2],paste(body(fi[[i]])))
        body(fi[[i]])<-gsub("p2",unlist(strsplit(modelList[[i]]$fct$name, "\\,|\\(|\\)"))[3],paste(body(fi[[i]])))
        f.all[[i]] <- fi[[i]]
      }
    }
    
    Body<-paste(mapply(FUN=function(x,y){paste(y, " * (",body(x),")")}, f.all, modelWeights), collapse = " + ")
    
    args <- as.character("dose")
    
    eval(parse(text = paste('g <- function(', args, ') { return(' , Body , "-", bmrScaled0, ')}', sep='')))
    
    if(identical(searchInterval,"dataBased")){
    LLimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]]))[2]/10000
    ULimit<-unique(sort(modelList[[1]]$data[[as.character(modelList[[1]]$call$formula)[[3]]]],decreasing=TRUE))[1]
    } else { 
      LLimit <- searchInterval[1]
      ULimit <- searchInterval[2]
    }
    BMD<-uniroot(g,interval=c(LLimit,ULimit))$root
    
    resMat<-matrix(c(BMD),1,1)
    colnames(resMat) <- c("BMD")
    rownames(resMat) <- c("")
  } else {
    curveNames <- colnames(modelList[[1]]$parmMat)
    uniqueCurveids <- unique(modelList[[1]]$dataList$curveid)
    nMods <- length(modelList)
    
    curveFcts <- lapply(modelList, 
                        function(model){ 
                          function(x) model$curve[[1]](x)
                        })
    
    getBMDi <- function(i){
      colIndex <- which(curveNames[i] == uniqueCurveids)
      
      g <- function(x){
        curveVal <- colSums(modelWeights * matrix(sapply(1:nMods, function(j) curveFcts[[j]](x)[, colIndex]), nrow = nMods, byrow = TRUE))
        curveVal - bmrScaled0[i]
      }
      
      if(identical(searchInterval,"dataBased")){
        LLimit <- unique(sort(modelList[[1]]$dataList$dose[modelList[[1]]$dataList$curveid == curveNames[i]]))[2]/10000
        ULimit <- unique(sort(modelList[[1]]$dataList$dose[modelList[[1]]$dataList$curveid == curveNames[i]],decreasing=TRUE))[1]
      } else { 
        LLimit <- searchInterval[i,1]
        ULimit <- searchInterval[i,2]
      }
      
      BMD <- uniroot(g,interval=c(LLimit,ULimit))$root
      BMD
    }
    
    BMD <- sapply(1:nCurves, getBMDi)
    
    resMat <- matrix(BMD, nrow = nCurves, ncol = 1)
    colnames(resMat) <- "BMD"
    rownames(resMat) <- colnames(modelList[[1]]$parmMat)
    g <- NULL
  }
  
  resBMD<-list(Results = resMat,
               MACurve = g)
  class(resBMD) <- "bmd"
  invisible(resBMD)
}


