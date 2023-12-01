qplotDrc <- function(x, add = FALSE, level = NULL, type = c("average", "all", "bars", "none", "obs", "confidence"), 
                      gridsize = 100, xtrans = "pseudo_log", xlab, xlim, 
                      ytrans = NULL, ylab, ylim, col = FALSE,
                      normal = FALSE, normRef = 1, confidence.level = 0.95){
  object <- x
  type <- match.arg(type)
  
  ## Determining logarithmic scales
  if ((xtrans == "log") || (xtrans == "pseudo_log"))
  {
    logX <- TRUE
  } else {
    logX <- FALSE
  }
  if(is.null(xtrans)){xtrans = "identity"}
  if(is.null(ytrans)){ytrans = "identity"}
  dataList <- object[["dataList"]]
  dose <- dataList[["dose"]]
  resp <- dataList[["origResp"]]
  curveid <- dataList[["curveid"]]
  plotid <- dataList[["plotid"]]
  
  
  if (identical(object[["type"]], "ssd")) 
  {
    dose <- unlist(with(dataList, tapply(dose, curveid, function(x){sort(x)}))[unique(dataList[["curveid"]])])
    resp <- unlist(with(dataList, tapply(dose, curveid, function(x){ppoints(x, 0.5)}))[unique(dataList[["curveid"]])])
  }
  
  ## Normalizing the response values
  if (normal)
  {
    names(resp) <- seq(length(resp))
    respList <- split(resp, curveid)
    
    respNorm <- mapply(normalizeLU, respList, 
                       as.list(as.data.frame(getLU(object)))[names(respList)], 
                       normRef = normRef, SIMPLIFY = F)
    
    resp <- do.call(c, unname(respNorm))[as.character(seq(length(resp)))]
  }
  
  if (!is.null(plotid))
  {       
    assayNoOld <- as.vector(plotid)  
  } else {
    assayNoOld <- as.vector(curveid)
  }
  uniAss <- unique(assayNoOld) 
  numAss <- length(uniAss) 

  
  doPlot <- is.null(level) || any(uniAss %in% level)
  if (!doPlot) {stop("Nothing to plot")}
  
  plotFct <- (object$"curve")[[1]]
  logDose <- (object$"curve")[[2]]
  naPlot <- ifelse(is.null(object$"curve"$"naPlot"), FALSE, TRUE)
 
  dlNames <- dataList[["names"]]
  doseName <- dlNames[["dName"]]
  respName <- dlNames[["orName"]]
  cNames <- dlNames[["cNames"]]
  
  if (missing(xlab)) {if (doseName == "") {xlab <- "Dose"} else {xlab <- doseName}}
  if (missing(ylab)) {if (respName == "") {ylab <- "Response"} else {ylab <- respName}}     
  
  ## Determining range of dose values
  if (missing(xlim)) 
  {
    xLimits <- c(min(dose), max(dose))
  } else {
    xLimits <- xlim  # if (abs(xLimits[1])<zeroEps) {xLimits[1] <- xLimits[1] + zeroEps}
  }
  
  if(logX){
    xLimits0 <- pmax(xLimits, 1e-8)
  }
  
  # Handling small dose values
  ## Constructing appropriate break on dose axis
  # if (!is.null(logDose))  # natural logarithm
  # {
  #   conLevel <- round(min(dose[is.finite(dose)])) - 1
  # } else {
  #   log10cl <- round(log10(min(dose[dose > 0]))) - 1
  #   conLevel <- 10^(log10cl)
  # }
  # 
  # if ((xLimits[1] < conLevel) && (logX || (!is.null(logDose))))
  # {
  #   xLimits[1] <- conLevel
  #   smallDoses <- (dose < conLevel)
  #   dose[smallDoses] <- conLevel
  # }
  # if (xLimits[1] >= xLimits[2]) {stop("Argument 'conLevel' is set too high")}

  ## Constructing dose values for plotting
  #    if (doseDim == 1) 
  #    {
  if ((is.null(logDose)) && (logX))
  {
    dosePts <- c(0,exp(seq(log(xLimits0[1]), log(xLimits0[2]), length = gridsize-1)))
    ## Avoiding that slight imprecision produces dose values outside the dose range
    ## (the model-robust predict method is sensitive to such deviations!)
    dosePts[1] <- max(xLimits[1],0)
    dosePts[gridsize] <- xLimits[2]           
  } else {
    dosePts <- seq(xLimits[1], xLimits[2], length = gridsize)
  }
  
  
  ## Finding minimum and maximum on response scale
  if (is.null(logDose)) 
  {
    plotMat <- plotFct(dosePts)
  } else {
    plotMat <- plotFct(logDose^(dosePts))
  }
  ## Normalizing the fitted values
  if (normal)
  {
    respList <- split(resp, curveid)
    plotMat <- mapply(normalizeLU, as.list(as.data.frame(plotMat)), 
                      as.list(as.data.frame(getLU(object))), 
                      normRef = normRef)            
    #        pmNew <- mapply(normalizeLU, as.list(as.data.frame(plotMat)), as.list(as.data.frame(getLU(object))))    
    #        print(pmNew)
    #        print(plotMat)
  }    
  #    numCol <- ncol(plotMat)
  
  maxR <- max(resp)
  options(warn = -1)  # suppressing warning in case maximum of NULL is taken
  maxPM <- apply(plotMat, 2, max, na.rm = TRUE)
  if (max(maxPM) > maxR) {maxPM <- maxPM[which.max(maxPM)]} else {maxPM <- maxR}
  options(warn=0)  
  
  if (missing(ylim)) 
  {
    yLimits <- NULL       
  } else {
    yLimits <- ylim
  }
  
  ## Cutting away original x values outside the limits
  eps1 <- 1e-8
  logVec <- !( (dose < xLimits[1] - eps1) | (dose > xLimits[2] + eps1) )
  dose <- dose[logVec]
  resp <- resp[logVec]
  #    assayNo <- assayNo[logVec]
  assayNoOld <- assayNoOld[logVec]
  
  ## Calculating predicted values for error bars
  if (identical(type, "bars"))
  {
    predictMat <- predict(object, interval = "confidence",
                          level = confidence.level)[, c("Lower", "Upper")]

    if(normal) {
      names(predictMat) <- seq(length(predictMat))
      predictList <- split(predictMat, curveid)
      predictMatListNorm <- mapply(normalizeLU, predictList,
                                   as.list(as.data.frame(getLU(object))),
                                   normRef = normRef,
                                   SIMPLIFY = F)
      predictMatNorm <- do.call(c, unname(predictMatListNorm))[as.character(seq(length(predictMat)))]
      predictMat<- matrix(predictMatNorm, ncol = 2)
    }
  }
  
  if (is.null(level)) 
  {
    level <- uniAss
  } else {
    level <- intersect(level, uniAss)
  }
  lenlev <- length(level)
  
  obsIndices <- which(assayNoOld %in% level)
  
  if((lenlev == 1) || (col == FALSE)) {
    colorAes <- "NULL"
  } else {
    colorAes <- "level"
  }
  if((lenlev == 1) || (col == TRUE)){
    shapeAes <- "NULL"
    linetypeAes <- "NULL"
  } else {
    shapeAes <- "level"
    linetypeAes <- "level"
  }
  
  obsLayer <- 
    switch(
      type, 
      "average" = ggplot2:::geom_point(aes(x = dose, y = resp, 
                                 col = eval(parse(text = colorAes)),
                                 shape = eval(parse(text = shapeAes))), 
                             data = data.frame(dose = dose[obsIndices], resp = resp[obsIndices], level = as.character(assayNoOld[obsIndices])) |> 
                                    group_by(dose,level) |> 
                                    tidyverse:::summarise(dose = first(dose),
                                              resp = mean(resp),
                                              level = first(level))),
      "bars"    = ggplot2:::geom_errorbar(aes(x = dose, ymin = Lower, ymax = Upper, 
                                    col = eval(parse(text = colorAes)),
                                    linetype = eval(parse(text = linetypeAes))),
                                data = data.frame(dose = dose[obsIndices], predictMat[obsIndices,], level = as.character(assayNoOld[obsIndices]))),
      "none"    = NULL, 
      "all"     = ggplot2:::geom_point(aes(x = dose, y = resp, 
                                 col = eval(parse(text = colorAes)),
                                 shape = eval(parse(text = shapeAes))),
                             data = data.frame(dose = dose[obsIndices], resp = resp[obsIndices], level = as.character(assayNoOld[obsIndices]))),
      "obs"     = ggplot2:::geom_point(aes(x = dose, y = resp, 
                                 col = eval(parse(text = colorAes)),
                                 shape = eval(parse(text = shapeAes))),
                             data = data.frame(dose = dose[obsIndices], resp = resp[obsIndices], level = as.character(assayNoOld[obsIndices])))
    )
  
  ## Plotting fitted curves
  curveLayer <- NULL
  if (!identical(type, "obs"))
  {
    curveLayer <- ggplot2:::geom_line(aes(x = x, y = y, 
                                col = eval(parse(text = colorAes)),
                                linetype = eval(parse(text = linetypeAes))), 
                            data = data.frame(x = rep(dosePts, lenlev),
                                              y = as.numeric(plotMat[,which(uniAss %in% level)]),
                                              level = rep(as.character(level), each = length(dosePts))))
  }
  
  # Confidence band
  if (identical(type, "confidence"))
  {
    confBandprLevel <- function(level0){
      newdata <- data.frame(DOSE=dosePts, CURVE=rep(level0, length(dosePts)))
      colnames(newdata) <- c(doseName, cNames)
      predictMat <- predict(object,
                            newdata=newdata,
                            interval = "confidence",
                            level=confidence.level)
      
      x <- c(dosePts, rev(dosePts))
      y <- c(predictMat[,"Upper"], rev(predictMat[,"Lower"]))
      
      confBandData <- data.frame(x,y, level = as.character(level0))
      ggplot2:::geom_polygon(aes(x,y, 
                       fill = eval(parse(text = colorAes)),
                       linetype = eval(parse(text = linetypeAes))), alpha = 0.5, data = confBandData)
    }
    confBandLayer <- lapply(level, confBandprLevel)
  } else {confBandLayer <- NULL}
  
  # Final plot
  if(!add){
    ggplot2:::ggplot() +
      confBandLayer +
      curveLayer +
      obsLayer +
      ggplot2:::scale_x_continuous(trans = xtrans, limits = xLimits) +
      ggplot2:::scale_y_continuous(trans = ytrans, limits = yLimits) +
      ggplot2:::labs(x = xlab, y = ylab, col = "", fill = "", shape = "", linetype = "")
  } else {
    list(
      confBandLayer = confBandLayer,
      curveLayer = curveLayer,
      obsLayer = obsLayer)
  }
}
