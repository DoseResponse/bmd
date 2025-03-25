\name{trendTest}
\alias{trendTest}
\title{Test for Trend in Relationship Between Two Variables}
\description{
  Conducts a test for trends between a numeric independent variable \code{x} and a numeric dependent variable \code{y} using specified statistical tests.
}
\usage{
trendTest(x, y, data, test = c("william", "shirley", "tukey"), level = 0.05)
}
\arguments{
  \item{x}{A numeric vector or the name of the independent variable (if \code{data} is provided).}
  \item{y}{A numeric vector or the name of the dependent variable (if \code{data} is provided).}
  \item{data}{An optional data frame containing the variables \code{x} and \code{y}. If provided, \code{x} and \code{y} should be column names in \code{data}.}
  \item{test}{A character string specifying the test to use. Must be one of \code{"william"}, \code{"shirley"}, or \code{"tukey"} (default).}
  \item{level}{Significance level for the test. Defaults to 0.05.}
}
\details{
  The function tests for a trend in the relationship between \code{x} and \code{y} based on the specified test:
  \itemize{
    \item \code{"william"}: Applies Williams' test to assess trend significance.
    \item \code{"shirley"}: Uses Shirley's test for trend analysis with ordered alternatives.
    \item \code{"tukey"}: Implements the Tukey trend test using multiple marginal models.
  }

  The direction of the trend (increasing or decreasing) is determined by the slope of the linear model \code{lm(y ~ x)}.
}
\value{
  A list with the following components:
  \item{p.values}{A numeric vector of p-values for the tests (if applicable).}
  \item{decisions}{A character vector indicating whether the trend is \code{"accept"} or \code{"reject"} based on the test results.}
  \item{acceptTrend}{A logical value indicating whether a trend is accepted (\code{TRUE}) or rejected (\code{FALSE}) based on the specified significance level.}
}
\examples{
# Example with custom data
x <- c(1, 2, 3, 4, 5)
y <- c(2, 4, 6, 8, 10)
result <- trendTest(x, y, test = "tukey")
print(result)

# Example with a data frame
data <- data.frame(x = c(1, 2, 3, 4, 5), y = c(10, 9, 8, 7, 6))
result <- trendTest("x", "y", data = data, test = "shirley")
print(result)
}
\seealso{
  \code{.williamsTest}, \code{.shirleyTest}, \code{.tukeytrendfit}
}
\references{
  Williams, D. A. (1971). "A test for differences between treatment means when several dose levels are compared with a zero dose control." Biometrics, 27(1), 103-117.
  Shirley, E. (1977). "A non-parametric equivalent of Williams' test for contrasting increasing dose levels of a treatment." Biometrics, 33(2), 386-389.
  Schaarschmidt, F. et al. (2021). "The Tukey trend test: Multiplicity adjustment using multiple marginal models" Biometrics, 78(2), 789-797.
}
\author{
  Jens Riis Baalkilde
}
\keyword{trend, statistical test}
