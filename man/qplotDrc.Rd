\name{qplotDrc}

\alias{qplotDrc}

\title{Plotting fitted dose-response curves using ggplot2}

\description{
  \code{qplotDrc} displays fitted curves and observations in the same plot window, 
  distinguishing between curves by different plot symbols and line types or colours using \code{ggplot2}.
}

\usage{
  
  qplotDrc(x, add = FALSE, level = NULL, type = c("average", "all", "bars", "none", "obs", "confidence"), 
            gridsize = 100, xtrans = "pseudo_log", xlab, xlim, 
            ytrans = NULL, ylab, ylim, col = FALSE,
            normal = FALSE, normRef = 1, confidence.level = 0.95)
}

\arguments{
  \item{x}{an object of class 'drc'.}
  \item{add}{logical. If TRUE then the functions returns a list of plot layers to be added to an already existing ggplot.} 
  \item{level}{vector of curve levels to plot. To plot only the curves specified by their names.}
  \item{type}{a character string specifying how to plot the data. There are currently
  5 options:  "average" (averages and fitted curve(s); default), 
  "none" (only the fitted curve(s)), "obs" (only the data points), 
  "all" (all data points and fitted curve(s)),
  "bars" (averages and fitted curve(s) with model-based standard errors (see Details)), and
  "confidence" (confidence bands for fitted curve(s)).}
  \item{gridsize}{numeric. Number of points in the grid used for plotting the fitted curves.}
  \item{xtrans}{Transformation to use on the x-axis. The default is "pesudo_log", which is a logarithmic 
  transformation with a smooth transition to 0. Other options include (among others) "log" and "identity". Can be overridden
  by adding a scale using the \code{scale_x_continuous} function.}
  \item{xlab}{an optional label for the x axis.}
  \item{xlim}{a numeric vector of length two, containing the lower and upper limit for the x axis.}
  \item{ytrans}{Transformation to use on the y-axis. The default is no transformation. Other options include (among others) "pseudo_log" and "log" and. Can be overridden
  by adding a scale using the \code{scale_y_continuous} function.}
  \item{ylab}{an optional label for the y axis.}
  \item{ylim}{a numeric vector of length two, containing the lower and upper limit for the y axis.}  
  \item{col}{logical. If TRUE default ggplot colours are used, can be overridden by \code{scale_color_manual}. 
  If FALSE (default) no colours are used.}
  \item{normal}{logical. If TRUE the plot of the normalized data and fitted curves are shown
  (for details see Weimer et al. (2012) for details).} 
  \item{normRef}{numeric specifying the reference for the normalization (default is 1).}  
  \item{confidence.level}{confidence level for error bars and confidence bands. Defaults to 0.95.}
}

\details{
  This function largely seeks to mimic the behaviour of the \code{plot} method for the \code{drc} package using 
  the \code{ggplot2} package for plotting.

  The use of \code{xlim} allows changing the range of the x axis, extrapolating the fitted dose-response curves. 
  Note that changing the range on the x axis may also entail a change of the range on the y axis.
  
  See \code{\link{colors}} for the available colours.

  Suitable labels are automatically provided.
  
  The model-based standard errors used for the error bars are calculated
  as the fitted value plus/minus the estimated error times the
  1-(alpha/2) quantile in the t distribution with degrees of freedom
  equal to the residual degrees of freedom for the model (or using a
  standard normal distribution in case of binomial and poisson data),
  where alpha=1-confidence.level. The standard errors are obtained using
  the predict method with the arguments interval = "confidence"
  and level=confidence.level.
}

\value{
  A \code{ggplot} object. If the option \code{add} is used, a list of \code{ggplot} layers is returned. 
}

\author{
  Jens Riis Baalkilde. Based on \code{plot.drc} by Christian Ritz and Jens C. Streibig with Contributions from Xiaoyan Wang and Greg Warnes. 
}

\references{
  Weimer, M., Jiang, X., Ponta, O., Stanzel, S., Freyberger, A., Kopp-Schneider, A. (2012)
  The impact of data transformations on concentration-response modeling.
  \emph{Toxicology Letters}, \bold{213}, 292--298.
}

%\note{ ~~further notes~~ }

%\seealso{\code{\link{plotraw}} plots the observations only.}

\examples{

## Fitting models to be plotted below
ryegrass.m1 <- drm(rootl~conc, data = ryegrass, fct = LL.4())
ryegrass.m2 <- drm(rootl~conc, data = ryegrass, fct = LL.3())  # lower limit fixed at 0

## Plotting observations and fitted curve for the first model
p <- qplotDrc(ryegrass.m1) 
p

## Add confidence region for the first model.
p + qplotDrc(ryegrass.m1, type="confidence", add = TRUE)$confBandLayer 

## Plot both models
p + qplotDrc(ryegrass.m2, add = TRUE)$curveLayer

## Fitting model to be plotted below
spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())

## Plot with no colours
qplotDrc(spinach.m1)

## Plot with default colours
qplotDrc(spinach.m1, col = TRUE)

## Plot with specified colours
qplotDrc(spinach.m1, col = TRUE) + 
  scale_color_manual(values = c("darkolivegreen4", "lightcoral", "goldenrod","darkslategray4", "plum4"))

## Plot of curves 1 and 2 only
qplotDrc(spinach.m1, level = c(1,4))

## Plot with confidence regions
qplotDrc(spinach.m1, col = TRUE, type = "confidence")

## Plot points and add confidence regions. Confidence regions are coloured by the "fill" aesthetic.
## Customising the x scale by adding a new scale.
qplotDrc(spinach.m1, col = TRUE, type = "confidence") +
  qplotDrc(spinach.m1, col = TRUE, type = "average", add = TRUE)$obsLayer + 
  scale_color_manual(values = c("darkolivegreen4", "lightcoral", "goldenrod","darkslategray4", "plum4"))+ 
  scale_fill_manual(values = c("darkolivegreen4", "lightcoral", "goldenrod","darkslategray4", "plum4")) +
  scale_x_continuous(trans = scales:::pseudo_log_trans(sigma = 0.2, base = exp(1)))


}
\keyword{ggplot}