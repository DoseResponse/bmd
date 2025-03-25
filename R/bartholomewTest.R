.bartholomewTest <- function (y, x, alternative = c("auto", "greater", "less"), number.of.bootstrap.samples = 1000, 
                               plot = NULL, seed = NULL) 
{
  alternative <- match.arg(alternative)
  if(alternative == "auto"){
    lm_slope <- lm(y ~ x)$coef[["x"]]
    slope <- ifelse(lm_slope > 0, "greater", "less")
  }
  
  # Summarised vector
  xFac <- factor(x)
  lm0 <- lm(y ~ xFac - 1)
  x <- y_mean <- summary(lm0)$coef[,"Estimate"]
  sigma <- summary(lm0)$coef[,"Std. Error"]
  
  # Start of LRT.trend function
  a <- 1/sigma^2
  xbar <- sum(a * x)/sum(a)
  k <- length(x)
  if(!is.null(seed)){
    set.seed(seed)
  }
  r <- matrix(stats::rnorm(k * number.of.bootstrap.samples), 
              ncol = k)
  r <- sweep(r, 2, sigma, "*")
  r <- sweep(r, 2, xbar, "+")
  
  if(alternative == "greater") {
    LRT.value.trend <- function (x, sigma) {
      a <- 1/sigma^2
      xbar <- sum(a * x)/sum(a)
      s <- seq_along(x)
      Atot <- cbind(s[-length(s)], s[-1])
      fit.ls1 <- isotone::activeSet(Atot, "LS", y = x, weights = a)
      LRT.increasing <- sum(a * (x - xbar)^2) - fit.ls1$fval
      return(LRT.increasing)
    }
  } else {
    LRT.value.trend <- function (x, sigma) {
      a <- 1/sigma^2
      xbar <- sum(a * x)/sum(a)
      s <- seq_along(x)
      Atot <- cbind(s[-1], s[-length(s)])
      fit.ls2 <- isotone::activeSet(Atot, "LS", y = x, weights = a)
      LRT.decreasing <- sum(a * (x - xbar)^2) - fit.ls2$fval
      return(LRT.decreasing) 
    }
  }
  
  L <- t(apply(r, 1, LRT.value.trend, sigma = sigma))
  obsL <- LRT.value.trend(x, sigma)
  
  STATISTIC = obsL
  PVAL <- mean(obsL <= L)
  
  RET <- list(statistic = STATISTIC, p.value = PVAL, alternative = alternative)
}
