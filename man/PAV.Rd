\name{PAV}
\alias{PAV}
\title{
Pool-adjacent-violators monotonizing
}
\description{The function monotonizes a sequence of probabilities or means based on the pool-adjacent-violators algorithm.  
}
\usage{
function(formula,data,type) 
}
\arguments{
  \item{formula}{an object of class "formula" expressing dose-respone relationship. Details of model specification are given under 'Details'}
  \item{data}{data frame containing the variables in the formula}
  \item{type}{character string specifying the type of data used in the model. "continuous" or "binomial"}
}

\details{
For details on how the pool-adjacent-violators algorithm is defined see Silvapulle and Sen (2004).

Formula should be specifed as in the form number/total ~ dose for binomial data and response ~ for continuous data.
}
\value{
A vector containing the monotonized sequence.
}
\references{
  Silvapulle, M. J. and Sen, P. K. (2004). Constrained statistical inference: order, inequality, and shape constraints. New York: John Wiley & Sons.
}

\author{ Signe M. Jensen
}
\note{
}

\seealso{
}

\keyword{nonparametric}
\keyword{isotonic regression}

\concept{Pool-adjacent-violators}