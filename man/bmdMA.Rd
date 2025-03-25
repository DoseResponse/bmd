\name{bmdMA}
\alias{bmdMA}
\title{
Model-averaged benchmark dose estimation
}
\description{Estimation of benchmark doses and benchmark dose lower limit based on model averaging from a user-defined list of dose response model fits 
}
\usage{
bmdMA(modelList, modelWeights, bmr, backgType = c("modelBased", "absolute", "hybridSD", "hybridPercentile"), backg = NA, 
def = c("excess", "additional", "relative", "extra", "added", "hybridExc", "hybridAdd", "point"),respTrans = c("none", "log", "sqrt"),
interval = c("delta","sandwich","inv", "profile"),
type = c("curve","bootstrap","Kang","Buckland"), bootstrapType, R=1000, bootInterval, level=0.95, stackingSeed, stackingSplits, display = TRUE, progressInfo = TRUE) 
}
\arguments{
  \item{modelList}{list of models of class \code{drc}}
  \item{modelWeights}{character string specifying the type of weights used, "AIC", "BIC" or "Stack", or a vector of the same length as the modelList with user defined weights}
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
    
  \item{backg}{numeric value specifying the background level. Defaults to 0 for "absolute" background risk for binomial response (1 for decreasing dose-response models), 2 SD for "hybridSD" background and 0.9 for "hybridpercentile"}
  
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
  
  \item{respTrans}{if the dose-response model is fitted with a transformed response, specifying this option ensures that the background level for the BMD is computed on the original scale. 
  Options include "none" (default), "log" (natural logarithm) and "sqrt"(square root).}
  
  \item{interval}{character string specifying the type of confidence interval to use for the individual models (no effect when type = "bootstrap" or "curve"): "delta" (default), "sandwich", "inv" or "profile"
  
  "delta" - BMDL is based on the lower limit of a Wald confidence interval based on the delta method
  
  "sandwich" - BMDL is based on the lower limit of a Wald confidence interval based on the delta method where the sandwich covariance matrix is used
  
  "inv" - BMDL is based on inverse regression, that is the dose associated with the upper limit of the point wise confidence interval of the dose-response curve. Not possible if type = "Buckland"
  
  "profile" - BMDL is based on the lower limit of a profile likelihood confidence interval. Not possible if type = "Buckland"}

  \item{type}{character string specifying how to estimate BMD and BMDL: "curve", "bootstrap", "Kang" or "Buckland" 
  
  "curve" - 
  
  "bootstrap" - 
  
  "Kang" - 
  
  "Buckland" - 
  }
  
  \item{bootstrapType}{character string indicating type of bootstrap sampling to be used if type="bootstrap". "nonparametric" (default), "semiparametric" or "parametric" (see details)}
  
  \item{bootInterval}{character string indicating how to estimate the bootstrap confidence intervals used to find BMDL type="bootstrap". "percentile" (default) or "BCa" (Bias corrected and adjusted)}
  
  \item{R}{number of bootstrap samples to use. Default is 1000}
  
  \item{level}{numeric value specifying the level of the confidence interval underlying BMDL. Default is 0.95. Note that the full CI is also computed. The full CI has level 1-2*(1-level)}
  
  \item{stackingSeed}{integer or NULL: Random seed to use in the data split in the estimation of the Stacking Weights when \code{modelWeights = "Stack"}. The global seed is reset to the initial value after estimation of the weights, so this option does not interfere with a globally set seed.}
  
  \item{stackingSplits}{integer or "LOO": When \code{modelWeights = "Stack"}, the Stacking Weights are estimated, which are based on V-fold cross-validation. The stackingSplits argument sets the number V of data splits used in the cross validation. The "LOO" (Leave one out) is a shortcut to setting V equal to the number of observations in the data set. }
  
  \item{display}{logical. If TRUE the results are displayed; otherwise they are not}
  
}

\details{
This package project is still under development. The aim to provide an R package 
  calculating the benchmark dose (BMD) and the lower limit of the corresponding 95\% confidence interval (BMDL)
  for continuous and quantal dose-response data for a range of dose-response model based on the available
  definitions of the benchmark dose concepts.
  
  Details on the implemented definitions and methods can be found in Crump (2002)
  
  Bootstrapping with the argument boot = "nonparametric" is done by sampling with replacement from the original data set. Bootstrapping with the argument boot = "parametric" is done by sampling from norm(mean(Y_i),sd(Y_0)), assuming equal variance between groups, in case of continuous data. For binomial data, each bootstrap data set is sampled from binom(N_i,Y_i/N_i). In case of Y_i = 0 or Y_i = N_i shrinkage is used to avoid that the resampling always produces 0 or 1, respectively. In this case data is sampled from binom(N_i,(Y_i+1/3)/(N_i+1/3)). Bootstrapping with argument boot = "semiparametric" is done by sampling with replacement from the residuals.

All sampling is made within dose groups.  
}
\value{
A list of three elements: Results contain the estimated BMD and BMDL, Boot.samples.used gives the number of boot samples that resulted in succesful estimations and were accordingly used  in the estimation of BMDL (and BMDU), Interval gives BMDL and BMDU, which is identical to the confidence interval for the percentile interval approach.
}
\references{
  Budtz-Jorgensen, E., Keiding, N., and Grandjean, P. (2001) 
  Benchmark Dose Calculation from Epidemiological Data,
  \emph{Biometrics} \bold{57}, 698--706.  

  Crump, K. (2002)
  Critical Issues in Benchmark Calculations from Continuous Data,
  \emph{Critical Reviews in Toxicology} \bold{32}, 133--153.
}

\author{ Signe M. Jensen and Jens Riis Baalkilde
}
\note{
}

\seealso{
}
\examples{
## Fitting 4 different two-parameter models to binomial data
deguelin.m1 <- drm(r/n~dose, weights=n, data=deguelin, fct=LL.2(), type="binomial")
deguelin.m2 <- drm(r/n~dose, weights=n, data=deguelin, fct=W1.2(), type="binomial")
deguelin.m3 <- drm(r/n~dose, weights=n, data=deguelin, fct=W2.2(), type="binomial")
deguelin.m4 <- drm(r/n~dose, weights=n, data=deguelin, fct=LN.2(), type="binomial")


## Model averaged BMD for 5% additional risk with estimated background risk and BMDL based on Buckland et al.
bmdMA(list(deguelin.m1,deguelin.m2,deguelin.m3,deguelin.m4), modelWeights="AIC", 0.05, 
      backgType = "modelBased", def = "additional",
      type = "Buckland")

## Model averaged BMD for 5% additional risk with estimated background risk and BMDL based on an average of the model curves
bmdMA(list(deguelin.m1,deguelin.m2,deguelin.m3,deguelin.m4), modelWeights="AIC", 0.05, 
      backgType = "modelBased", def = "additional",
      type = "curve", bootstrapType = "parametric", bootInterval = "percentile", R=100)


## Fitting 4 different two-parameter models to binomial data
ryegrass.m1<-drm(rootl~conc, data=ryegrass, fct=LL.4())
ryegrass.m2<-drm(rootl~conc, data=ryegrass, fct=W1.4())
ryegrass.m3<-drm(rootl~conc, data=ryegrass, fct=W2.4())
ryegrass.m4<-drm(rootl~conc, data=ryegrass, fct=LN.4())

## Model-averaged BMD and bootstrap BMDL for bmr=5% and using the hybrid approach to estimate the background risk.  
bmdMA(list(ryegrass.m1,ryegrass.m2,ryegrass.m3,ryegrass.m4), modelWeights="AIC", bmr=0.05, 
      backgType = "hybridSD",
      def = "hybridAdd", type = "bootstrap",bootstrapType = "nonparametric", bootInterval = "percentile")


## Model-averaged BMD using the Stacking Weights
bmdMA(list(ryegrass.m1,ryegrass.m2,ryegrass.m3,ryegrass.m4), modelWeights="Stack", bmr=0.05, 
      backgType = "hybridSD",
      def = "hybridAdd", type = "bootstrap",bootstrapType = "nonparametric", bootInterval = "percentile")

}
\keyword{model averaging}
\keyword{nonlinear}
\keyword{bootstrap}

\concept{BMD BMDL benchmark dose-response model averaging}