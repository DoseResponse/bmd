\name{bmdOrdinal}
\alias{bmdOrdinal}
\title{
Benchmark dose estimation for ordinal dose-response models
}
\description{Estimation of benchmark doses and benchmark dose lower limit from ordinal dose response model fits 
}
\usage{
bmdOrdinal(object, bmr, backgType = c("modelBased", "absolute"), backg = NA,
           def = c("excess", "additional", "point"), 
           interval = c("delta", "sandwich", "profile", "bootstrap"), level = 0.95, 
           R = 500, bootType = c("nonparametric", "parametric", "model", "hierarchical"), 
           display = TRUE, progressInfo = TRUE) 
}
\arguments{
  \item{object}{object of class \code{drcOrdinal}}
  \item{bmr}{numeric value of benchmark response level for which to calculate the benchmark dose}
  \item{backgType}{character string specifying how the background level is specified. The options are "modelBased" and "absolute".

  "modelBased" - the background level is obtained from the model as the level for dose 0: 
  p0 = f(0) 
  
  "absolute" - the background level is specified by the user through the backg argument.}
    
  \item{backg}{numeric value specifying the background level. Defaults to 0 for "absolute" background risk.}
  
  \item{def}{character string specifying the definition of the benchmark dose to use in the calculations. "excess", "additional" and "point" are available for ordinal response.
  
  "excess" - BMR is defined as: BMR = (f(BMD) - p0)/(1 - p0). 
  Works for binomial response. BMR should be between 0 and 1.  
  
  "additional" - BMR is defined as: BMR = f(BMD) - p0.
  Works for binomial response. BMR should be between 0 and 1.  
  
  "point" - The response level for which to find BMD is directly defined through the BMR level: BMR = f(BMD). Works for binomial, count and continuous response data.}

  \item{interval}{character string specifying the type of confidence interval to use: "delta" (default), "sandwich", "profile" or "bootstrap" 
  
  "delta" - BMDL is based on the lower limit of a Wald confidence interval based on the delta method
  
  "sandwich" - BMDL is based on the lower limit of a Wald confidence interval based on the delta method, where the sandwich covariance matrix is used
  
  "bootstrap" - BMDL is based on bmd estimates on resampled data sets. The type of bootstrap applied can be specified by the "bootType" argument. Default is "nonparametric".
  
  "none" - No confidence interval is computed.}
  
  \item{level}{numeric value specifying the levle of the confidence interval underlying BMDL. Default is 0.95}
  
  \item{R}{integer specifying the number of data sets resampled from the original data set when computing the confidence interval by bootstrap.}
  
  \item{bootType}{character string specifying the resampling procedure for the data sets used for bootstrap.
  
  "nonparametric" - Bootstrapping is done by sampling with replacement from the original data set.
  
  "parametric" - Bootstrapping is done by sampling from a multinomial distribution with probabilites given by the number of observations for each level. If all observations for one group are in the same level, shrinkage is used to avoid that the resampling always produces that particular level. In this case data is sampled from a multinomial distribution with probabilities (1/|K|^2)/(N_i + 1/|K|) for the levels with 0 observations, and (N_i + 1/|K|^2)/(N_i + 1/|K|), where N_i is the number of observatiions for the particular group, and |K| is the number of levels.
  
  "model" - New data sets are resampled from the fitted model.
  
  "hierarchical" - If the experiment design is hierarchical (i.e. observations are grouped according to a categorical variable, for instance litter), a special resampling procedure is suggested (Blessinger 2020)
  
  All resampling is done within the dose values. 
  }
  
  \item{display}{logical. If TRUE the results are displayed; otherwise they are not}
  
  \item{progressInfo}{logical. If TRUE, a progressbar is displayed when calculating bootstrap confidence intervals}
}

\details{
The aim to provide an R package 
  calculating the benchmark dose (BMD) and the lower limit of the corresponding 95\% confidence interval (BMDL)
  for continuous, quantal, quantal and ordinal dose-response data for a range of dose-response models based on the available
  definitions of the benchmark dose concepts.
  
  Details on the implemented definitions and methods can be found in Crump (2002)
}
\value{
A list of four elements: Results contain the estimated BMD and BMDL, interval gives the lower (BMDL) and upper (BMDU) end of the confidence interval of BMD. 
}
\references{
  Budtz-Jorgensen, E., Keiding, N., and Grandjean, P. (2001) 
  Benchmark Dose Calculation from Epidemiological Data,
  \emph{Biometrics} \bold{57}, 698--706.  

  Crump, K. (2002)
  Critical Issues in Benchmark Calculations from Continuous Data,
  \emph{Critical Reviews in Toxicology} \bold{32}, 133--153.
  
  Blessinger, T. D., Euling, S. Y., Wang, L., Hogan, K. A., Cai, C., Klinefelter, G., Saillenfait, A.M. (2020)
  \emph{Environment International} \bold{134}
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
data(guthion)

guthionS <- subset(guthion, trt == "S")

guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), 
                          weights = "total", dose = "dose", data = guthionS, fct = LL.2())
bmdOrdinal(guthionS.LL, bmr = 0.1, backgType = "modelBased", def = "excess", R = 50)
}
\keyword{models}
\keyword{nonlinear}

\concept{BMD BMDL benchmark dose-response}