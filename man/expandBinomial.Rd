\name{expandBinomial}
\alias{expandBinomial}
\title{
Utility function
}
\description{Utility function for expanding binomial data
}
\usage{
expandBinomial(data, number, total, dose, curveid = character(0))
}

\arguments{
  \item{data}{a data.frame}
  \item{number}{the name of the column in the data set containing the number of affected individuals per dose level}
  \item{total}{the name of the column in the data set containing the total number of individuals per dose level}
  \item{dose}{the name of the column in the data set containing the dose values}
  \item{curveid}{the name of the column in the data set specifying the curveid (if available)}
}

\details{
The aim to provide an R package 
  calculating the benchmark dose (BMD) and the lower limit of the corresponding 95\% confidence interval (BMDL)
  for continuous and quantal dose-response data for a range of dose-response models based on the available
  definitions of the benchmark dose concepts.

}
\value{
data.frame
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
}
\keyword{models}
\keyword{nonlinear}

\concept{BMD BMDL benchmark dose-response}