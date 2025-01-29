\name{drmHetVar}
\alias{drmHetVar}
\title{
Dose response modeling with heterogeneous variance
}
\description{Add a heterogeneous variance structure to an existing \code{drc} object.
}
\usage{
drmHetVar(object, var.formula)
}

\arguments{
  \item{object}{dose-response model of class \code{drcHetVar}}
  \item{var.formula}{one-sided formula specifying the dependance of the dose values and/or the fitted values on the point-wise standard error}
}

\details{
The aim to provide an R package 
  calculating the benchmark dose (BMD) and the lower limit of the corresponding 95\% confidence interval (BMDL)
  for continuous and quantal dose-response data for a range of dose-response models based on the available
  definitions of the benchmark dose concepts.
  
  REFERENCES TO BE ADDED/WRITTEN
}
\value{
dose-response model with a heterogeneous variance structure of class \code{drcHetVar}. 

The primary objective is to use this model for benchmark dose estimation based on the hybrid method with a heterogeneous variance structure.

A plot method is available, which can be useful for assessing the fit of the variance structure.
}
\references{
}

\author{ Signe M. Jensen and Jens Riis Baalkilde
}
\note{
}

\seealso{
}
\examples{
library(drc)
library(drcData)
library(bmd)
# install.packages("gridExtra") # OPTIONAL - USED FOR PLOTTING A drcHetVar OBJECT.

# ryegrass data
ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
set.seed(123)
ryegrass.LL.4.hetVar <- drmHetVar(ryegrass.LL.4, ~ fitted + I(fitted^2))
plot(ryegrass.LL.4.hetVar)
bmdHetVar(ryegrass.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE)

# barley data
barley.LL.4 <- drm(weight ~ Dose, data = barley, fct = LL.4())
set.seed(123)
barley.LL.4.hetVar <- drmHetVar(barley.LL.4, ~ fitted + I(fitted^2))
plot(barley.LL.4.hetVar)
bmdHetVar(barley.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE)
bmdHetVar(barley.LL.4.hetVar, bmr = 0.1, backgType = "hybridSD", backg = 1, def = "hybridExc", R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE)

# GiantKelp data
GiantKelp.LL.4 <- drm(tubeLength ~ dose, data = GiantKelp, fct = LL.4())
set.seed(123)
GiantKelp.LL.4.hetVarSq <- drmHetVar(GiantKelp.LL.4, ~ fitted + I(fitted^2))
plot(GiantKelp.LL.4.hetVarSq)
bmdHetVar(GiantKelp.LL.4.hetVarSq, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE)
bmdHetVar(GiantKelp.LL.4.hetVarSq, bmr = 0.1, backgType = "hybridSD", backg = 1, def = "hybridExc", R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE)

GiantKelp.LL.4.hetVarLogSq <- drmHetVar(GiantKelp.LL.4, ~ log(dose+1) + I(log(dose+1)^2))
plot(GiantKelp.LL.4.hetVarLogSq)
bmdHetVar(GiantKelp.LL.4.hetVarLogSq, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE)
bmdHetVar(GiantKelp.LL.4.hetVarLogSq, bmr = 0.1, backgType = "hybridSD", backg = 1, def = "hybridExc", R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE)



}
\keyword{models}
\keyword{nonlinear}

\concept{BMD BMDL benchmark dose-response}