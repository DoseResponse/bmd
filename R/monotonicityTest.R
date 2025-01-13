monotonicityTest <- function(x, y, data, test = c("jonckheere", "bartholomew"), level = 0.05, ...){ # , "drc", "quad"
  if(!missing(data)){
    x <- data[[x]]
    y <- data[[y]]
  }
  
  
  xFac <- factor(x)
  
  lm_alternative <- lm(y ~ x)$coef[["x"]]
  alternative <- ifelse(lm_alternative > 0, "greater", "less")
  
  test <- match.arg(test)
  if(test == "jonckheere"){
    p.value <- .jonckheereTest(x = y, g = x, alternative = alternative)$p.value
    names(p.value) <- NULL
    acceptMonotonicity = p.value < level
  }
  
  if(test == "bartholomew"){
    p.value <- .bartholomewTest(y = y, x = x, alternative = alternative, ...)$p.value
    names(p.value) <- NULL
    acceptMonotonicity = p.value < level
  }
  
  # if(test == "drc"){
  #   capture.output({
  #     p.value <- .drcMonotonicityTest(y = y, x = x, alternative = alternative, ...)$p.value
  #   })
  #   names(p.value) <- NULL
  #   acceptMonotonicity = p.value > level
  # }
  # 
  # if(test == "quad"){
  #   p.value <- .quadMonotonicityTest(y = y, x = x, ...)$p.value
  #   names(p.value) <- NULL
  #   acceptMonotonicity = p.value > level
  # }
  
  list(p.value = p.value, acceptMonotonicity = acceptMonotonicity)
}
