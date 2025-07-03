#' Help funtion to bmdBoot
#' 
#' 
#' 
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


