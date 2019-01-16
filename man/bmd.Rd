\name{bmd}
\alias{bmd}
\title{
Benchmark dose estimation
}
\description{Estimation of benchmark doses and benchmark dose lower limit from dose response model fits 
}
\usage{
bmd(object, bmr, backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"), backg, controlSD,
def = c("excess", "additional", "relative", "extra", "added", "hybridExc", "hybridAdd", "point"), 
              interval = c("delta","sandwich","inv"), display = FALSE) 
}
\arguments{
  \item{object}{object of class \code{drc}}
  \item{bmr}{numeric value of benchmark response level for which to calculate the benchmark dose}
  \item{backgType}{character string specifying how the background level is specified. For binomial data the options are "modelBased" and "absolute". For continuous data the options are "absolute", "hybridSD" and "hybridPercentile"

  "modelBased" - the background risk is obtained from the model as the risk for dose 0: 
  p0 = f(0) 
  
  "absolute" - the background risk is specified by the user through the backg argument:
  p0 = backg for binomial response and for the "relative", "extra" and "added" definition for continuous response. 
  p0 = 1 - phi((back - f(0))/sigma) for "hybridExc" and "hybridAdd" definitions.
  
  "hybridSD" - the background risk is specified by the user in terms of number of SDs from the mean of the control group.
  p0 = 1 - phi(((backg*sigma + f(0)) - f(0))/sigma) = 1 - phi(backg),
  where phi is the normal distribution function and sigma is the SD for the control group. 
  
  "hybridPercentile" - the background risk is specified by the user in terms of percentile from the control group distribution (assuming a normal distribution).
  p0 = 1 - phi((x0 - f(0))/sigma) = 1 - backg.
  where x0 is the level for which the response is considered adverse, phi is the normal distribution function and sigma is the SD for the control group}
    
  \item{backg}{numeric value specifying the background level. Defaults to 0 for "absolute" background risk for binomial response (1 for decreasing dose-response models), 2 SD for "hybridSD" background and 0.9 for "hybridPercentile"}
  
  \item{controlSD}{numeric value specifying the standard deviation of the control group (used in the hybrid approach). If not specified the SD for the control group will be estimated as the square root of the residual variance (assuming variance homogeneity).}
  
  \item{def}{character string specifying the definition of the benchmark dose to use in the calculations. "excess" , "additional" and "point" are for binomial response whereas "relative", "extra", "added", "hybridExc" (excess hybrid), "hybridAdd" (additional hybrid), and "point" are for continuous response
  
  "excess" - BMR is defined as: BMR = (f(BMD) - p0)/(1 - p0). 
  Works for binomial response. BMR should be between 0 and 1.  
  
  "additional" - BMR is defined as: BMR = f(BMD) - p0.
  Works for binomial response. BMR should be between 0 and 1.  
  
  "point" - The response level for which to find BMD is directly defined through the BMR level: BMR = f(BMD). Works for both binomial and continuous response 
  
  "relative" - BMR is defined as: BMR = (f(BMD) - p0)/p0. 
  Works for continuous response 
  
  "extra" - BMR is defined as: BMR = (f(BMD) - p0)/(f(Inf) - p0).
  Works for continuous response
  
  "added" - BMR is defined as: BMR= f(BMD) + p0.
  Works for continuous response
  
  "hybridExc" - BMR is defined as: BMR = (1 - phi((x0 - f(BMD))/sigma) - p0)/ (1- p0), 
  where x0 is the level for which the response is considered adverse, phi is the normal distribution function and sigma is the SD for the control group. 
  Works for continuous response
  
  "hybridAdd" - BMR is defined as: BMR = 1 - phi((x0 - f(BMD))/sigma) - p0, 
  where x0 is the level for which the response is considered adverse, phi is the normal distribution function and sigma is the SD for the control group. 
  Works for continuous response}
  
  \item{interval}{character string specifying the type of confidence interval to use: "delta" (default), "sandwich" or "inv" 
  
  "delta" - BMDL is based on the lower limit of a Wald confidence interval based on the delta method
  
  "sandwich" - BMDL is based on the lower limit of a Wald confidence intervals based on the delta method with robust standard errors from the sandwich covariance matrix 
  
  "inv" - BMDL is based on inverse regression, that is the dose associated with the upper limit of the pointwise confidence interval of the dose-response curve}
  
  \item{display}{logical. If TRUE the results are displayed; otherwise they are not}
}

\details{
This package project is still under development. The aim to provide an R package 
  calculating the benchmark dose (BMD) and the lower limit of the corresponding 95\% confidence interval (BMDL)
  for continuous and quantal dose-response data for a range of dose-response model based on the available
  definitions of the benchmark dose concepts.
  
  Details on the implemented definitions and methods can be found in Crump (2002)
}
\value{
A matrix with two columns, one containing BMD and the other containing BMDL.
}
\references{
  Budtz-Jorgensen, E., Keiding, N., and Grandjean, P. (2001) 
  Benchmark Dose Calculation from Epidemiological Data,
  \emph{Biometrics} \bold{57}, 698--706.  

  Crump, K. (2002)
  Critical Issues in Benchmark Calculations from Continuous Data,
  \emph{Critical Reviews in Toxicology} \bold{32}, 133--153.
}

\author{ Signe M. Jensen
}
\note{
}

\seealso{
}
\examples{
## Fitting log-logistic two-parameter model to binomial data
deguelin.m1 <- drm(r/n~dose, weights=n, data=deguelin, fct=LL.2(), type="binomial")

## BMD for 5% additional risk with estimated background risk
bmd(deguelin.m1, 0.05, backgType = "modelBased", def = "additional")

## BMD for 10% additional risk with 2% background risk
bmd(deguelin.m1, 0.1, backg = 0.02 , backgType = "absolute", def = "additional")

## BMD for 5% excess risk and background 0
bmd(deguelin.m1, 0.05, backg = 0, backgType = "absolute", def = "excess")

## Dose resulting in 12% risk
bmd(deguelin.m1, 0.12, def = "point")

## Benchmark doses for a continuous response
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## BMD as the dose resulting in a 5% change relative to the mean background level
bmd(ryegrass.m1, 0.05, backgType = "modelBased", def = "relative", display = TRUE)

## BMD using the hybrid method, background risk is 2 SD, hybrid definition using excess risk
bmd(ryegrass.m1, 0.05, backg = 2, backgType = "hybridSD", def = "hybridAdd", display = TRUE)

}
\keyword{models}
\keyword{nonlinear}

\concept{BMD BMDL benchmark dose-response}