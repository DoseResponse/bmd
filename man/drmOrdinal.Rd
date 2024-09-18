\name{drmOrdinal}
\alias{drmOrdinal}
\title{
Fitting ordinal dose-response models
}
\description{A model fitting function for analysis of ordinal dose-response data.
}
\usage{
drmOrdinal(levels, dose, weights, data, fct, p.epsilon = 10^(-16))
}
\arguments{
  \item{levels}{a character vector of the levels of the response variable in increasing order, as they appear in the supplied data set}
  
  \item{dose}{a character string specifying the column with the dose values in the supplied data set}
  
  \item{weights}{a character string specifying the column containing the number of observations pr. group}
  
  \item{blocks}{a character string specifying the column containing the blocks of the experiment, if available (optional)}
  
  \item{data}{a dataframe containing the observations}
  
  \item{fct}{a list with three or more elements specifying the non-linear
  function, the accompanying self starter function, the names of the parameter in the non-linear function and, 
  optionally, the first and second derivatives as well, for the individual fitted curves. For more information see the help page for the "drm" function in the "drc" package.
  Currently available functions for ordinal dose-response models include, among others, the log-logistic models 
  \code{\link{LL.4}}, the log-normal model \code{\link{LN.4}} and the two Weibull models \code{\link{W1.4}} and \code{\link{W2.4}}. Use
  \code{\link{getMeanFunctions}} for a full list.}

  \item{p.epsilon}{numeric value specifying the lower bound for the probabilites for each level returned by the function pFun created when fitting the model. Default value is 10^(-16)}
}

\details{
 This functions fits a dose-response model for ordinal dose-response data by fitting a series of binomial dose-response models. 
}
\value{
  An object of (S3) class 'drcOrdinal'. 
}
\references{
  # ADD REFERENCES
}

\author{ Signe M. Jensen and Jens Riis Baalkilde
}
\note{
}

\seealso{
}
\examples{
library(drcData)
data(guthion)
guthionS <- subset(guthion, trt == "S")

guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = LL.2())

plot(guthionS.LL, xlim = c(15,55)) # uses ggplot

}
\keyword{models}
\keyword{nonlinear}

\concept{BMD BMDL benchmark dose-response}