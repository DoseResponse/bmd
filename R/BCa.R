#' Calculate BCa (Bias-Corrected and Accelerated) Bootstrap Confidence Interval
#'
#' @description
#' This function calculates the BCa confidence interval for a given statistic using bootstrap and jackknife estimates.
#' 
#' @details
#' It computes the BCa confidence interval for a statistic using bootstrap and jackknife estimates.
#' The BCa method adjusts for both bias and skewness in the bootstrap distribution and generally
#' provides better coverage than standard bootstrap confidence intervals.
#'
#' @param obs numeric; The observed value of the statistic computed from the original sample
#' @param data data.frame or matrix; The original dataset from which bootstrap samples are drawn
#' @param bootSample numeric vector; Bootstrap replicates of the statistic. Must be of length R 
#'        where R is the number of bootstrap replicates
#' @param bootjack numeric vector; Jackknife replicates of the statistic. Must be of length n 
#'        where n is the number of observations in data
#' @param level numeric; Confidence level between 0 and 1 (e.g., 0.95 for 95% confidence interval)
#'
#' @return A numeric vector of length 2 containing:
#' \itemize{
#'   \item lower: The lower confidence limit
#'   \item upper: The upper confidence limit
#' }
#'
#'
#' @examples \dontrun{
#' # Example with mean as the statistic
#' data <- rnorm(100)
#' obs_mean <- mean(data)
#' 
#' # Generate bootstrap samples
#' R <- 1000
#' boot_means <- replicate(R, mean(sample(data, replace = TRUE)))
#' 
#' # Generate jackknife samples
#' jack_means <- sapply(1:length(data), 
#'                      function(i) mean(data[-i]))
#' 
#' # Calculate 95% CI
#' BCa(obs_mean, data, boot_means, jack_means, 0.95)
#' }
#'
#' @references
#' Efron, B. (1987). Better Bootstrap Confidence Intervals. 
#' _Journal of the American Statistical Association_, 82(397), 171-185.
#'
#' DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
#' _Statistical Science_, 11(3), 189-212.
#'
#' @export
BCa <- function(obs, data, bootSample, bootjack, level){
  R <- length(bootSample)
  b <- qnorm((sum(bootSample > obs)+sum(bootSample==obs)/2)/R)
  
  n <- nrow(data) 
  n1 <- n-1 
  obsn <- obs*n
  pv <- i <- 0 
  while(i < n){
    i = i+1 
    pv[i] = obsn-n1*bootjack[i]
  }
  pv<-pv[!is.na(pv)]
  
  je <- mean(pv)-pv
  a <- sum(je^3)/(6*sum(je^2))^(3/2)
  
  alpha <- (1-level)*2 
  z <- qnorm(c(alpha/2,1-alpha/2)) # Std. norm. limits
  p <- pnorm((z-b)/(1-a*(z-b))-b) # correct & convert to proportions
  
  quantile(bootSample,p=p) # ABC percentile lims.      
}









