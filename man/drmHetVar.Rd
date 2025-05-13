\name{drmHetVar}
\alias{drmHetVar}
\title{
Dose response modeling with heterogeneous variance
}
\description{Fit a dose-response model with heterogeneous variance dependending on dose level.
}
\usage{
drmHetVar(formula, var.formula, data, fct)
}

\arguments{
  \item{formula}{formula for the dose-response relationship}
  \item{var.formula}{one-sided formula specifying the dependance of the dose values and/or the fitted values on the point-wise standard error}
  \item{data}{data.frame containing the observations}
  \item{fct}{a list with three or more elements specifying the non-linear function, the accompanying self starter function, the names of the parameter in the non-linear function and, optionally, the first and second derivatives as well as information used for calculation of ED values. Currently available functions include, among others, the four- and five-parameter log-logistic models LL.4, LL.5 and the Weibull model W1.4. Use drc::getMeanFunctions for a full list.}
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

The primary objective is to use this model for benchmark dose estimation based on the hybrid method with a heterogeneous variance structure. This can be done using the \code{bmdHetVar} function.

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
set.seed(123)
ryegrass.LL.4.hetVar <- drmHetVar(rootl ~ conc, ~ fitted + I(fitted^2), data = ryegrass, fct = LL.4())
plot(ryegrass.LL.4.hetVar)
bmdHetVar(ryegrass.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 500, level = 0.95, progressInfo = TRUE, display = TRUE)
bmdHetVar(ryegrass.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 500, level = 0.95, 
          bootType = "parametric", progressInfo = TRUE, display = TRUE) # parametric bootstrap

# barley data
set.seed(123)
barley.LL.4.hetVar <- drmHetVar(weight ~ Dose, ~ fitted + I(fitted^2), data = barley, fct = LL.4())
plot(barley.LL.4.hetVar)
bmdHetVar(barley.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 500, level = 0.95, progressInfo = TRUE, display = TRUE)
bmdHetVar(barley.LL.4.hetVar, bmr = 0.1, backgType = "hybridSD", backg = 1, def = "hybridExc", R = 500, level = 0.95, progressInfo = TRUE, display = TRUE)

# GiantKelp data
set.seed(123)
GiantKelp.LL.4.hetVarSq <- drmHetVar(tubeLength ~ dose, ~ fitted + I(fitted^2), data = GiantKelp, fct = LL.4())
plot(GiantKelp.LL.4.hetVarSq)
bmdHetVar(GiantKelp.LL.4.hetVarSq, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 500, level = 0.95, progressInfo = TRUE, display = TRUE)
bmdHetVar(GiantKelp.LL.4.hetVarSq, bmr = 0.1, backgType = "hybridSD", backg = 1, def = "hybridExc", R = 500, level = 0.95, progressInfo = TRUE, display = TRUE)

GiantKelp.LL.4.hetVarLogSq <- drmHetVar(tubeLength ~ dose, ~ log(dose+1) + I(log(dose+1)^2), data = GiantKelp, fct = LL.4())
plot(GiantKelp.LL.4.hetVarLogSq)
bmdHetVar(GiantKelp.LL.4.hetVarLogSq, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 500, level = 0.95, progressInfo = TRUE, display = TRUE)
bmdHetVar(GiantKelp.LL.4.hetVarLogSq, bmr = 0.1, backgType = "hybridSD", backg = 1, def = "hybridExc", R = 500, level = 0.95, progressInfo = TRUE, display = TRUE)



}
\keyword{models}
\keyword{nonlinear}

\concept{BMD BMDL benchmark dose-response}