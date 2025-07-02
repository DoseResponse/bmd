\name{plot.bmd}

\alias{plot}

\title{Plotting benchmark doses}

\description{
  \code{plot.bmd} displays benchmark dose values with options to plot confidence intervals as well.
}

\usage{
  plot(x, ..., interval = c("BMDL", "twosided", "none"))
}

\arguments{
  \item{x}{an object of class 'bmd'.}
  \item{...}{arguments to be passed on to \code{plot.drc}, if \code{add = FALSE}}
  \item{interval}{option to plot only the lower limit of the confidence interval for the benchmark dose ("BMDL", default), both limit
  of the confidence interval ("twosided"), or no confidence interval ("none").}
}

\details{
  This function is a simple function to plot benchmark dose values along with the fitted curve.
}

\value{
  Creates a plot. No value returned.
}

\author{
  Jens Riis Baalkilde.
}

%\references{
%}

%\note{ ~~further notes~~ }

%\seealso{\code{\link{plotraw}} plots the observations only.}

\examples{
library(drc)
library(drcData)

## Fitting model and calculating BMD. 
model <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
bmd0 <- bmd(model, bmr = 0.1, backgType = "modelBased", def = "relative")

## Plotting
plot(bmd0)

# Plot both limits of confidence interval
plot(bmd0, interval = "twosided")

# Pass argument to plot.bmd to plot confidence band around curve
plot(bmd0, type = "confidence")


}
\keyword{ggplot}