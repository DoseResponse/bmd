\name{qplotBmd}

\alias{qplotBmd}

\title{Plotting benchmark doses using ggplot2}

\description{
  \code{qplotBmd} displays benchmark dose values with options to plot confidence intervals as well using \code{ggplot2}.
}

\usage{
  qplotBmd(x, ..., interval = c("BMDL", "twosided", "none"), col = FALSE, add = FALSE)
}

\arguments{
  \item{x}{an object of class 'bmd'.}
  \item{...}{arguments to be passed on to qplotDrc, if \code{add = FALSE}}
  \item{interval}{option to plot only the lower limit of the confidence interval for the benchmark dose ("BMDL", default), both limit
  of the confidence interval ("twosided"), or no confidence interval ("none").}
  \item{col}{logical. If TRUE then multiple curves specified by "curveid" in the dose-response model are distinguised by colours rather than point shapes and line types} 
  \item{add}{logical. If TRUE then the functions returns a list of plot layers to be added to an already existing ggplot.} 
}

\details{
  This function is a simple function to plot benchmark dose values 
}

\value{
  A \code{ggplot} object. If the option \code{add} is used, a list of \code{ggplot} layers is returned. 
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

# Plot
qplotBmd(bmd0, interval = "twosided", add = FALSE)

qplotDrc(model,type="confidence") +
  qplotBmd(bmd0, interval = "twosided", add = TRUE)

qplotBmd(bmd0, interval = "twosided", add = FALSE)


}
\keyword{ggplot}