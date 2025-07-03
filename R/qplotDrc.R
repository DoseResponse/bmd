#' Plotting fitted dose-response curves using ggplot2
#' 
#' \code{qplotDrc} displays fitted curves and observations in the same plot
#' window, distinguishing between curves by different plot symbols and line
#' types or colours using \code{ggplot2}.
#' 
#' This function largely seeks to mimic the behaviour of the \code{plot} method
#' for the \code{drc} package using the \code{ggplot2} package for plotting.
#' 
#' The use of \code{xlim} allows changing the range of the x axis,
#' extrapolating the fitted dose-response curves.  Note that changing the range
#' on the x axis may also entail a change of the range on the y axis.
#' 
#' See \code{\link{colors}} for the available colours.
#' 
#' Suitable labels are automatically provided.
#' 
#' The model-based standard errors used for the error bars are calculated as
#' the fitted value plus/minus the estimated error times the 1-(alpha/2)
#' quantile in the t distribution with degrees of freedom equal to the residual
#' degrees of freedom for the model (or using a standard normal distribution in
#' case of binomial and poisson data), where alpha=1-confidence.level. The
#' standard errors are obtained using the predict method with the arguments
#' interval = "confidence" and level=confidence.level.
#' 
#' @param x an object of class 'drc'.
#' @param add logical. If TRUE then the functions returns a list of plot layers
#' to be added to an already existing ggplot.
#' @param level vector of curve levels to plot. To plot only the curves
#' specified by their names.
#' @param type a character string specifying how to plot the data. There are
#' currently 5 options: "average" (averages and fitted curve(s); default),
#' "none" (only the fitted curve(s)), "obs" (only the data points), "all" (all
#' data points and fitted curve(s)), "bars" (averages and fitted curve(s) with
#' model-based standard errors (see Details)), and "confidence" (confidence
#' bands for fitted curve(s)).
#' @param gridsize numeric. Number of points in the grid used for plotting the
#' fitted curves.
#' @param xtrans Transformation to use on the x-axis. The default is
#' "pesudo_log", which is a logarithmic transformation with a smooth transition
#' to 0. Other options include (among others) "log" and "identity". Can be
#' overridden by adding a scale using the \code{scale_x_continuous} function.
#' @param xlab an optional label for the x axis.
#' @param xlim a numeric vector of length two, containing the lower and upper
#' limit for the x axis.
#' @param ytrans Transformation to use on the y-axis. The default is no
#' transformation. Other options include (among others) "pseudo_log" and "log"
#' and. Can be overridden by adding a scale using the \code{scale_y_continuous}
#' function.
#' @param ylab an optional label for the y axis.
#' @param ylim a numeric vector of length two, containing the lower and upper
#' limit for the y axis.
#' @param col logical. If TRUE default ggplot colours are used, can be
#' overridden by \code{scale_color_manual}.  If FALSE (default) no colours are
#' used.
#' @param normal logical. If TRUE the plot of the normalized data and fitted
#' curves are shown (for details see Weimer et al. (2012) for details).
#' @param normRef numeric specifying the reference for the normalization
#' (default is 1).
#' @param confidence.level confidence level for error bars and confidence
#' bands. Defaults to 0.95.
#' @return A \code{ggplot} object. If the option \code{add} is used, a list of
#' \code{ggplot} layers is returned.
#' @author Jens Riis Baalkilde. Based on \code{plot.drc} by Christian Ritz and
#' Jens C. Streibig with Contributions from Xiaoyan Wang and Greg Warnes.
#' @references Weimer, M., Jiang, X., Ponta, O., Stanzel, S., Freyberger, A.,
#' Kopp-Schneider, A. (2012) The impact of data transformations on
#' concentration-response modeling.  \emph{Toxicology Letters}, \bold{213},
#' 292--298.
#' @keywords ggplot
#' @examples
#' 
#' library(drc)
#' library(drcData)
#' library(ggplot2)
#' 
#' ## Fitting models to be plotted below
#' ryegrass.m1 <- drm(rootl~conc, data = ryegrass, fct = LL.4())
#' ryegrass.m2 <- drm(rootl~conc, data = ryegrass, fct = LL.3())  # lower limit fixed at 0
#' 
#' ## Plotting observations and fitted curve for the first model
#' p <- qplotDrc(ryegrass.m1) 
#' p
#' 
#' ## Add confidence region for the first model.
#' p + qplotDrc(ryegrass.m1, type="confidence", add = TRUE)$confBandLayer 
#' 
#' ## Plot both models
#' p + qplotDrc(ryegrass.m2, add = TRUE)$curveLayer
#' 
#' ## Fitting model to be plotted below
#' spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())
#' 
#' ## Plot with no colours
#' qplotDrc(spinach.m1)
#' 
#' ## Plot with default colours
#' qplotDrc(spinach.m1, col = TRUE)
#' 
#' ## Plot with specified colours
#' qplotDrc(spinach.m1, col = TRUE) + 
#'   scale_color_manual(values = c("darkolivegreen4", "lightcoral", "goldenrod",
#'                                 "darkslategray4", "plum4"))
#' 
#' ## Plot of curves 1 and 2 only
#' qplotDrc(spinach.m1, level = c(1,4))
#' 
#' ## Plot with confidence regions
#' qplotDrc(spinach.m1, col = TRUE, type = "confidence")
#' 
#' ## Plot points and add confidence regions. Confidence regions are coloured by the "fill" aesthetic.
#' ## Customising the x scale by adding a new scale.
#' qplotDrc(spinach.m1, col = TRUE, type = "confidence") +
#'   qplotDrc(spinach.m1, col = TRUE, type = "average", add = TRUE)$obsLayer + 
#'   scale_color_manual(values = c("darkolivegreen4", "lightcoral", "goldenrod",
#'                                 "darkslategray4", "plum4"))+ 
#'   scale_fill_manual(values = c("darkolivegreen4", "lightcoral", "goldenrod",
#'                                "darkslategray4", "plum4")) +
#'   scale_x_continuous(trans = scales:::pseudo_log_trans(sigma = 0.2, base = exp(1)))
#' 
#' 
#' @export
qplotDrc <- function(x, add = FALSE, level = NULL, type = c("average", "all", "bars", "none", "obs", "confidence"), 
                      gridsize = 250, xtrans = "pseudo_log", xlab, xlim, 
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
  
  # Constructing dose values for plotting
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
    if(is.null(object$objList)){
      predictMat <- predict(object, interval = "confidence",
                            level = confidence.level)[, c("Lower", "Upper")]
    } else {
      predictMatList <- lapply(object$objList, function(x){
        predict(x, interval = "confidence", level = confidence.level)[, c("Lower", "Upper")]
      })
      predictMat <- do.call(rbind, predictMatList)
    }

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
      "average" = geom_point(aes(x = dose, y = resp, 
                                 col = eval(parse(text = colorAes)),
                                 shape = eval(parse(text = shapeAes))), 
                             data = data.frame(dose = dose[obsIndices], resp = resp[obsIndices], level = as.character(assayNoOld[obsIndices])) |> 
                                    group_by(dose,level) |> 
                                    summarise(dose = first(dose),
                                              resp = mean(resp),
                                              level = first(level))),
      "bars"    = geom_errorbar(aes(x = dose, ymin = Lower, ymax = Upper, 
                                    col = eval(parse(text = colorAes)),
                                    linetype = eval(parse(text = linetypeAes))),
                                data = data.frame(dose = dose[obsIndices], predictMat[obsIndices,], level = as.character(assayNoOld[obsIndices]))),
      "none"    = NULL, 
      "all"     = geom_point(aes(x = dose, y = resp, 
                                 col = eval(parse(text = colorAes)),
                                 shape = eval(parse(text = shapeAes))),
                             data = data.frame(dose = dose[obsIndices], resp = resp[obsIndices], level = as.character(assayNoOld[obsIndices]))),
      "obs"     = geom_point(aes(x = dose, y = resp, 
                                 col = eval(parse(text = colorAes)),
                                 shape = eval(parse(text = shapeAes))),
                             data = data.frame(dose = dose[obsIndices], resp = resp[obsIndices], level = as.character(assayNoOld[obsIndices])))
    )
  
  ## Plotting fitted curves
  curveLayer <- NULL
  if (!identical(type, "obs"))
  {
    curveLayer <- geom_line(aes(x = x, y = y, 
                                col = eval(parse(text = colorAes)),
                                linetype = eval(parse(text = linetypeAes))), 
                            data = data.frame(x = rep(dosePts, lenlev),
                                              y = as.numeric(plotMat[,which(uniAss %in% level)]),
                                              level = rep(as.character(level), each = length(dosePts))))
  }
  
  # Confidence band
  if (identical(type, "confidence"))
  {
    if(is.null(object$objList)){
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
        geom_polygon(aes(x,y, 
                         fill = eval(parse(text = colorAes)),
                         linetype = eval(parse(text = linetypeAes))), alpha = 0.5, data = confBandData)
      }
    }
    else {
      confBandprLevel <- function(level0){
        newdata <- data.frame(DOSE=dosePts)
        colnames(newdata) <- c(doseName)
        predictMat <- predict(object$objList[[level0]],
                              newdata=newdata,
                              interval = "confidence",
                              level=confidence.level)
        
        x <- c(dosePts, rev(dosePts))
        y <- c(predictMat[,"Upper"], rev(predictMat[,"Lower"]))
        
        confBandData <- data.frame(x,y, level = as.character(level0))
        geom_polygon(aes(x,y, 
                         fill = eval(parse(text = colorAes)),
                         linetype = eval(parse(text = linetypeAes))), alpha = 0.5, data = confBandData)
      }
    }
    confBandLayer <- lapply(level, confBandprLevel)
  } else {confBandLayer <- NULL}
  
  # Final plot
  if(!add){
    ggplot() +
      confBandLayer +
      curveLayer +
      obsLayer +
      scale_x_continuous(trans = xtrans, limits = xLimits) +
      scale_y_continuous(trans = ytrans, limits = yLimits) +
      labs(x = xlab, y = ylab, col = "", fill = "", shape = "", linetype = "")
  } else {
    list(
      confBandLayer = confBandLayer,
      curveLayer = curveLayer,
      obsLayer = obsLayer)
  }
}
