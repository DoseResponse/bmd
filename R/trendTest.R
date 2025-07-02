trendTest <- function(x, y, data, test = c("william", "shirley", "tukey"), level = 0.05){
  if(!missing(data)){
    x <- data[[x]]
    y <- data[[y]]
  }
  
  lm_slope <- lm(y ~ x)$coef[["x"]]
  slope <- ifelse(lm_slope > 0, "greater", "less")
  
  test <- match.arg(test)
  if(test == "william"){
    res <- .williamsTest(y, x, alternative = slope)
    p.values <- NULL
    decisions <- ifelse(res$statistic > res$crit.value, "accept", "reject")
    acceptTrend <- any(res$statistic > res$crit.value)
  }
  
  if(test == "shirley"){
    res <- .shirleyTest(y, x, alternative = slope, method = "look-up")
    p.values <- NULL
    decisions <- ifelse(res$statistic > res$crit.value, "accept", "reject")
    acceptTrend <- any(res$statistic > res$crit.value)
  }
  
  if(test == "tukey"){
    if(!requireNamespace("multcomp")){
      stop('package "multcomp" must be installed to use tukey trend test')
    }
    fitw <- lm(y ~ x)
    ttw <- .tukeytrendfit(y, x)
    res <- summary(multcomp::glht(model=ttw$mmm, linfct=ttw$mlf))
    
    p.values <- as.numeric(res$test$pvalues)
    names(p.values) <- names(res$test$tstat)
    decisions <- ifelse(p.values < level, "accept", "reject")
    acceptTrend <- any(p.values < level)
  }
  
  list(p.values = p.values, decisions = decisions, acceptTrend = acceptTrend)
}