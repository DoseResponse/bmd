\name{bmdBoot}
\alias{bmdBoot}
\title{
Benchmark dose estimation using bootstrap 
}
\description{The function estimates BMD and BMDL using bootstrap based on parametric dose-response models.  
}
\usage{
bmdBoot(object, bmr, R = 1000, bootType = "nonparametric", bmdType = "orig",
        backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"), 
        backg = NA, controlSD = NA,
        def = c("excess", "additional", "relative", "extra",
                "added", "hybridExc", "hybridAdd", "point"), 
        respTrans = c("none", "log", "sqrt"), bootInterval = c("percentile", "BCa"),
        display = TRUE, level = 0.95)
}
\arguments{
  \item{object}{object of class \code{drc}}
  \item{bmr}{numeric value of benchmark response level for which to calculate the benchmark dose}
  \item{R}{number of bootstrap samples. Default is 1000}
  \item{bootType}{character string specifying type of bootstrap used. "nonparametric" (default), "semiparametric" or "parametric". "Semiparametric is only available for continuous data and "nonparametric" is at present the only option for count data. See details below}
  \item{bmdType}{Type of estimate for BMD. Default is "orig" the bmd estimate from the original data set. Other choices are "mean" - the mean of the bootstrap estimates, or "median" - the median of the bootstrap estimates}
  \item{backgType}{character string specifying how the background level is specified. For binomial and count data the options are "modelBased" and "absolute". For continuous data the options are "modelBased", "absolute", "hybridSD" and "hybridPercentile"

  "modelBased" - the background level is obtained from the model as the level for dose 0: 
  p0 = f(0) 
  
  "absolute" - the background level is specified by the user through the backg argument:
  p0 = backg for binomial response and for the "relative", "extra" and "added" definition for continuous response. 
  p0 = 1 - phi((back - f(0))/sigma) for "hybridExc" and "hybridAdd" definitions.
  
  "hybridSD" - the background risk is specified by the user in terms of number of SDs from the mean of the control group.
  p0 = 1 - phi(((backg*sigma + f(0)) - f(0))/sigma) = 1 - phi(backg),
  where phi is the normal distribution function and sigma is the SD for the control group. 
  "hybridPercentile" - the background risk is specified by the user in terms of percentile from the control group distribution (assuming a normal distribution).
  p0 = 1 - phi((x0 - f(0))/sigma) = 1 - backg.
  where x0 is the level for which the response is considered adverse, phi is the normal distribution function and sigma is the SD for the control group}
    
  \item{backg}{numeric value specifying the background level. Defaults to 0 for "absolute" background level for binomial response (1 for decreasing dose-response models), 2 SD for "hybridSD" background and 0.9 for "hybridpercentile"}
  
  \item{controlSD}{numeric value specifying the standard deviation of the control group (used in the hybrid approach). If not specified the SD for the control group will be estimated as the square root of the residual variance (assuming variance homogeneity).}
  
  \item{def}{character string specifying the definition of the benchmark dose to use in the calculations. "excess" , "additional" and "point" are for binomial response whereas "relative", "extra", "added", "hybridExc" (excess hybrid), "hybridAdd" (additional hybrid), and "point" are for continuous response. "relative", "extra", and "point" are for count response data.
  
  "excess" - BMR is defined as: BMR = (f(BMD) - p0)/(1 - p0). 
  Works for binomial response. BMR should be between 0 and 1.  
  
  "additional" - BMR is defined as: BMR = f(BMD) - p0.
  Works for binomial response. BMR should be between 0 and 1.  
  
  "point" - The response level for which to find BMD is directly defined through the BMR level: BMR = f(BMD). Works for binomial, count and continuous response 
  
  "relative" - BMR is defined as: BMR = (f(BMD) - p0)/p0. 
  Works for count and continuous response 
  
  "extra" - BMR is defined as: BMR = (f(BMD) - p0)/(f(Inf) - p0).
  Works for count and continuous response
  
  "added" - BMR is defined as: BMR= f(BMD) + p0.
  Works for count and continuous response
  
  "hybridExc" - BMR is defined as: BMR = (1 - phi((x0 - f(BMD))/sigma) - p0)/ (1- p0), 
  where x0 is the level for which the response is considered adverse, phi is the normal distribution function and sigma is the SD for the control group. 
  Works for continuous response
  
  "hybridAdd" - BMR is defined as: BMR = 1 - phi((x0 - f(BMD))/sigma) - p0, 
  where x0 is the level for which the response is considered adverse, phi is the normal distribution function and sigma is the SD for the control group. 
  Works for continuous response}
  
  \item{respTrans}{if the dose-response model is fitted with a transformed response, specifying this option ensures that the background level for the BMD is computed on the original scale. 
  Options include "none" (default), "log" (natural logarithm) and "sqrt"(square root).}
  
  \item{bootInterval}{character string indicating type of bootstrap interval used for estimating BMDL. "Percentile" (default) or "BCa" (Bias corrected and adjusted)}

  \item{display}{logical. If TRUE the results are displayed; otherwise they are not}
  
  \item{level}{numeric value specifying the levle of the confidence interval underlying BMDL. Default is 0.95}

}

\details{
BMDL is defined as the 5th percentile in the bootstrap distribution.

Bootstrapping with the argument boot = "nonparametric" is done by sampling with replacement from the original data set. Bootstrapping with the argument boot = "parametric" is done by sampling from norm(mean(Y_i),sd(Y_0)), assuming equal variance between groups, in case of continuous data. For binomial data, each bootstrap data set is sampled from binom(N_i,Y_i/N_i). In case of Y_i = 0 or Y_i = N_i shrinkage is used to avoid that the resampling always produces 0 or 1, respectively. In this case data is sampled from binom(N_i,(Y_i+1/4)/(N_i+1/2)). Bootstrapping with argument boot = "semiparametric" is done by sampling with replacement from the residuals.

All sampling is made within dose groups.  

}
\value{
A list of three elements: Results contain the estimated BMD and BMDL, bootEst is a vector af all of the estimated BMD values from each bootstrap sample, Interval gives BMDL and BMDU, which is identical to the confidence interval for the percentile interval approach. 
}
\references{

}

\author{ Signe M. Jensen
}
\note{
}

\seealso{
}
\examples{
## Data on root length in ryegrass after exposure to ferulic acid
require(drc)
require(drcData)
data(ryegrass)

ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## BMD using the hybrid method, background risk is 2 SD, hybrid definition using excess risk
bmdBoot(ryegrass.m1, 0.05, backgType = "hybridSD", def = "hybridAdd", R = 50)

## BMD from the same definitions but using parametric bootstrap
bmdBoot(ryegrass.m1, 0.05, backgType = "hybridSD", def = "hybridAdd", bootType="parametric",R = 50)

}
\keyword{bootstrap}

\concept{BMD BMDL benchmark dose-response}