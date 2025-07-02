\name{monotonicityTest}
\alias{monotonicityTest}
\title{Test for Monotonicity of Relationship Between Two Variables}
\description{
  Conducts a test for monotonicity between a numeric independent variable \code{x} and a numeric dependent variable \code{y} using specified statistical tests.
}
\usage{
monotonicityTest(x, y, data, test = c("jonckheere", "bartholomew"), level = 0.05, ...)
}
\arguments{
  \item{x}{A numeric vector or the name of the independent variable (if \code{data} is provided).}
  \item{y}{A numeric vector or the name of the dependent variable (if \code{data} is provided).}
  \item{data}{An optional data frame containing the variables \code{x} and \code{y}. If provided, \code{x} and \code{y} should be column names in \code{data}.}
  \item{test}{A character string specifying the test to use. Must be one of \code{"jonckheere"} (default) or \code{"bartholomew"}.}
  \item{level}{Significance level for the test. Defaults to 0.05.}
  \item{...}{Additional arguments passed to the underlying test functions.}
}
\details{
  The function tests the monotonicity of the relationship between \code{x} and \code{y} based on the specified test:
  \itemize{
    \item \code{"jonckheere"}: Uses the Jonckheere-Terpstra test to assess monotonic trends.
    \item \code{"bartholomew"}: Uses Bartholomew's test to assess monotonicity.
  }

  The direction of the monotonicity (increasing or decreasing) is determined by the sign of the coefficient from a simple linear model \code{lm(y ~ x)}.
}
\value{
  A list with the following components:
  \item{p.value}{The p-value of the test.}
  \item{acceptMonotonicity}{A logical value indicating whether monotonicity is accepted (\code{TRUE}) or rejected (\code{FALSE}) based on the specified significance level.}
}
\examples{
# Example with custom data
x <- c(1, 2, 3, 4, 5)
y <- c(2, 4, 6, 8, 10)
result <- monotonicityTest(x, y, test = "jonckheere")
print(result)

data <- data.frame(x = c(1, 2, 3, 4, 5), y = c(10, 9, 8, 7, 6))
result <- monotonicityTest("conc", "rootl", data = drcData::ryegrass, test = "bartholomew")
print(result)
}

\author{Jens Riis Baalkilde}

\seealso{
  \code{.jonckheereTest}, \code{.bartholomewTest}
}

\references{
  A. R. Jonckheere (1954). "A Distribution-Free k-Sample Test Against Ordered Alternatives."
  D. J. Bartholomew (1961). "Ordered tests in the analysis of variance."
  OECD (2006). Rapport No. 54, Annexes.
}
\keyword{monotonicity, trend test}
