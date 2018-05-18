\name{bmdBoot}
\alias{bmdBoot}
\title{
Benchmark dose estimation using bootstrap 
}
\description{The function estimates BMD and BMDL using bootstrap based on parametric dose-response models.  
}
\usage{
bmd(object, bmr, R, boot, bmdType, backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"),
    backg=NA, def = c("excess", "additional", "relative", "extra", "added", "hybridExc", "hybridAdd", "point")) 
}
\arguments{
  \item{object}{object of class \code{drc}}
  \item{bmr}{numeric value of benchmark response level for which to calculate the benchmark dose}
  \item{R}{number of bootstrap samples. default is 1000}
  \item{boot}{character string specifying type of bootstrap used. "resample" (default) or "pseudorandom". See details below}
  \item{bmdType}{Type of estimate for BMD. Default is "orig" the bmd estimate from the original data set. Other choises are "mean" - the mean of the bootstrap estiamtes, or "median" - the median of the bootstrap estiamtes}
  \item{backgType}{character string specifying how the background level is specified. For binomial data the options are "modelBased" and "absolute". For continuous data the options are "absolute", "hybridSD" and "hybridPercentile"

  "modelBased" - the background risk is obtained from the model as the risk for dose 0: 
  p0 = f(0) 
  
  "absolute" - the background risk is specified by the user through the backg argument:
  p0 = backg for binomial response and for the "relative", "extra" and "added" definiton for continuous response. 
  p0 = 1 - phi((back - f(0))/sigma) for "hybridExc" and "hybridAdd" definitions.
  
  "hybridSD" - the background risk is specified by the user in terms of number of SDs from the mean of the control group.
  p0 = 1 - phi(((backg*sigma + f(0)) - f(0))/sigma) = 1 - phi(backg),
  where phi is the normal distribution function and sigma is the SD for the control group. 
  "hybridPercentile" - the background risk is specified by the user in terms of percentile from the control group distribution (assuming a normal distribution).
  p0 = 1 - phi((x0 - f(0))/sigma) = 1 - backg.
  where x0 is the level for which the response is considered adverse, phi is the normal distribution function and sigma is the SD for the control group}
    
  \item{backg}{numeric value specifying the background level. Defaults to 0 for "absolute" background risk for binomail response (1 for decreasing dose-response models), 2 SD for "hybridSD" background and 0.9 for "hybridpercentile"}
  
  \item{def}{character string specifiying the definition of the benchmark dose to use in the calculations. "excess" , "additional" and "point" are for binomial response whereas "relative", "extra", "added", "hybridExc" (excess hybrid), "hybridAdd" (additional hybrid), and "point" are for continuous response
  
  "excess" - BMR is defined as: BMR = (f(BMD) - p0)/(1 - p0). 
  Works for binomial response. BMR should be between 0 and 1.  
  
  "additional" - BMR is defined as: BMR = f(BMD) - p0.
  Works for binomial response. BMR should be between 0 and 1.  
  
  "point" - The response level for which to find BMD is directly defined through the BMR level: BMR = f(BMD). Works for both binomial and continious response 
  
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

}

\details{
BMDL is defined as the 5th percentile in the bootstrap distribution.

Bootstrapping with the argument boot = "resample" is done by sampling with replacement from the original data set. Bootstrapping with the argument boot = "pseudorandom" is done by sampling from norm(mean(Y_i),sd(Y_0)), assuming equal variance between groups, in case of continuous data. For binomial data, each bootstrap data set is sampled from binom(N_i,Y_i/N_i). In case of Y_i = 0 or Y_i = N_i shinkage is used to avoid that the resampling always produces 0 or 1, respectively. In this case data is sampled from binom(N_i,(Y_i+1/3)/(N_i+1/3)). 

All sampling is made within dose groups.  

}
\value{
A matrix with two columns, one containing BMD and the other containing BMDL.
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
data(ryegrass)

ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## BMD using the hybrid method, background risk is 2 SD, hybrid definition using excess risk
bmdBoot(ryegrass.m1, 0.05, backgType = "hybridSD", def = "hybridAdd")

## BMD from the same definitions but using parametric bootstrap
bmdBoot(ryegrass.m1, 0.05, backgType = "hybridSD", def = "hybridAdd", boot="pseudorandom")

}
\keyword{bootstrap}

\concept{BMD BMDL benchmark dose-response}