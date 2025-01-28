\name{bmd}
\alias{bmd}
\title{
Benchmark dose estimation with heterogeneous variance
}
\description{Estimation of benchmark doses and benchmark dose lower limit based on the hybrid method from dose response model fits with the option to specify a heterogeneous variance structure, where the variance depends on the dose level and/or the fitted values
}
\usage{
bmdHetVar(object, bmr, backgType = c("absolute", "hybridSD", "hybridPercentile"), 
          backg = NA, def = c("hybridExc", "hybridAdd"), interval = c("boot", "none"),
          R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE)
}

\arguments{
  \item{object}{dose-response model with a heterogeneous variance structure of class \code{drcHetVar}}
  \item{bmr}{numeric value of benchmark response level for which to calculate the benchmark dose}
  \item{backgType}{character string specifying how the background level is specified. For binomial data the options are "modelBased" and "absolute". For continuous data the options are "modelBased","absolute", "hybridSD" and "hybridPercentile". For count data (Poisson, negbin1 or negbin2) the options are "modelBased" and "absolute".

  "absolute" - the background level is specified by the user through the backg argument:
  p0 = 1 - phi((back - f(0))/sigma(0)) for "hybridExc" and "hybridAdd" definitions.
  
  "hybridSD" - the background risk is specified by the user in terms of number of SDs from the mean of the control group.
  p0 = 1 - phi(((backg*sigma(0) + f(0)) - f(0))/sigma(0)) = 1 - phi(backg),
  where phi is the normal distribution function and sigma(0) is the SD for the control group. 
  
  "hybridPercentile" - the background risk is specified by the user in terms of percentile from the control group distribution (assuming a normal distribution).
  p0 = 1 - phi((x0 - f(0))/sigma(0)) = 1 - backg.
  where x0 is the level for which the response is considered adverse, phi is the normal distribution function and sigma(0) is the SD for the control group}
    
  \item{backg}{numeric value specifying the background level. Defaults to 0 for "absolute" background risk for binomial response (1 for decreasing dose-response models), 2 SD for "hybridSD" background and 0.9 for "hybridPercentile"}
  
  \item{def}{character string specifying the definition of the benchmark dose to use in the calculations. "hybridExc" (excess hybrid), "hybridAdd" (additional hybrid), available.
  
  "hybridExc" - BMR is defined as: BMR = (1 - phi((x0 - f(BMD))/sigma(BMD)) - p0)/ (1- p0), 
  where x0 is the level for which the response is considered adverse, phi is the normal distribution function and sigma(BMD) is the SD at the benchmark dose. 
  
  "hybridAdd" - BMR is defined as: BMR = 1 - phi((x0 - f(BMD))/sigma(BMD)) - p0, 
  where x0 is the level for which the response is considered adverse, phi is the normal distribution function and sigma(BMD) is the SD at the benchmark dose. }
  
  \item{interval}{character string specifying the type of confidence interval to use: "boot" (default) or "none"
  
  "boot" - BMDL is based on non-parametric bootstrapping. 
  
  "none" - no confidence interval is computed.} 
  
  \item{R}{Number of bootstrap samples. Ignored if \code{interval = "none"}}
  
  \item{level}{numeric value specifying the levle of the confidence interval underlying BMDL. Default is 0.95}
  
  \item{progressInfo}{logical. If TRUE, progress info is be printed while bootstrap confidence intervals are estimated. Default is TRUE.}
  
  \item{display}{logical. If TRUE the results are displayed; otherwise they are not}
}

\details{
The aim to provide an R package 
  calculating the benchmark dose (BMD) and the lower limit of the corresponding 95\% confidence interval (BMDL)
  for continuous and quantal dose-response data for a range of dose-response models based on the available
  definitions of the benchmark dose concepts.
  
  REFERENCES TO BE ADDED/WRITTEN
}
\value{
A list of five elements: Results contain the estimated BMD and BMDL, bmrScaled is the response value corresponding to the BMD, interval gives the lower (BMDL) and upper (BMDU) end of the confidence interval of BMD, sigmaFun is the estimated standard deviation function, and model returns the supplied model.
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
# install.packages("gridExtra") # OPTIONAL - USED FOR PLOTTING sigmaFun

# ryegrass data
ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
set.seed(123)
ryegrass.LL.4.hetVar <- drmHetVar(ryegrass.LL.4, ~ fitted + I(fitted^2))
plot(ryegrass.LL.4.hetVar)
bmdHetVar(ryegrass.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE)
bmdHetVar(ryegrass.LL.4.hetVar, bmr = 0.1, backgType = "hybridSD", backg = 1, def = "hybridExc", R = 1000, level = 0.95, progressInfo = TRUE, display = TRUE)

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