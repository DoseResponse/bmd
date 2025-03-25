.shirleyTest <- function (x, g, alternative = c("auto", "greater", "less"), nperm = 10000, ...) 
{
  # if (is.list(x)) {
  #   if (length(x) < 2L) 
  #     stop("'x' must be a list with at least 2 elements")
  #   DNAME <- deparse(substitute(x))
  #   x <- lapply(x, function(u) u <- u[complete.cases(u)])
  #   k <- length(x)
  #   l <- sapply(x, "length")
  #   if (any(l == 0)) 
  #     stop("all groups must contain data")
  #   g <- factor(rep(1:k, l))
  #   alternative <- x$alternative
  #   nperm <- x$nperm
  #   method <- x$method
  #   x <- unlist(x)
  # }
  # else {
  if (length(x) != length(g)) {
    stop("'x' and 'g' must have the same length")
  }
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  OK <- complete.cases(x, g)
  x <- x[OK]
  g <- g[OK]
  if (!all(is.finite(g))) {
    stop("all group levels must be finite")
  }
  g.old <- g
  g <- factor(g)
  k <- nlevels(g)
  if (k < 2) {
    stop("all observations are in the same group")
  }
  # }
  alternative <- match.arg(alternative)
  if(alternative == "auto"){
    lm0 <- lm(x ~ g.old)$coef[["g.old"]]
    alternative <- ifelse(lm0 > 0, "greater", "less")
  }
  
  xold <- x
  if (alternative == "less") {
    x <- -x
  }
  nj <- tapply(x, g, length)
  k <- nlevels(g)
  kk <- k - 1
  if (kk > 10) 
    stop("Critical t-values are only available for up to 10 dose levels.")
  N <- sum(nj)
  compfn <- function(x, ix, g, nj) {
    k <- length(nj)
    ti <- rep(NA, k)
    x <- x[ix]
    for (i in 2:k) {
      N <- sum(nj[1:i])
      r <- rank(x[1:N])
      gg <- g[1:N]
      Rj <- tapply(r, gg, mean)
      t <- table(r)
      names(t) <- NULL
      T <- sum((t^3 - t)/(12 * (N - 1)))
      Vi <- N * (N + 1)/12 - T
      u <- 2:i
      j <- u
      enum <- sapply(j, function(j) sum(nj[j:i] * Rj[j:i]))
      denom <- sapply(j, function(j) sum(nj[j:i]))
      ti[i] <- (max(enum/denom) - Rj[1])/sqrt(Vi * (1/nj[i] + 
                                                      1/nj[1]))
    }
    return(ti[2:k])
  }
  l <- 1:N
  STATISTIC <- compfn(x, l, g, nj)
  
  extrapolFN <- function(Tki, beta, r, c) {
    out <- Tki - 0.01 * beta * (1 - r/c)
    return(out)
  }
  df <- 1000000
  c <- nj[1]
  r <- nj[2:k]
  nrows <- nrow(williams.tk005) # PMCMRplus:::TabCrit$williams.tk005
  Tkdf <- numeric(kk)
  dft <- as.numeric(williams.tk005$rowname) # PMCMRplus:::TabCrit$williams.tk005 # as.numeric(rownames(williams.tk005))
  xx <- c(2:6, 8, 10)
  for (i in 2:kk) {
    if (i <= 6 | i == 8 | i == 10) {
      yt <- williams.tk005[, paste0("X", i)] # PMCMRplus:::TabCrit$williams.tk005
      yb <- williams.beta005[, paste0("X", i)] # PMCMRplus:::TabCrit$williams.beta005
    }
    else {
      yt <- sapply(1:nrows, function(j) {
        approx(x = xx, y = williams.tk005[j,-1], xout = i)$y # PMCMRplus:::TabCrit$williams.tk005
      })
      yb <- sapply(1:nrows, function(j) {
        approx(x = xx, y = williams.beta005[j,-1], xout = i)$y # PMCMRplus:::TabCrit$williams.tk005
      })
    }
    tt <- approx(x = dft, y = yt, xout = df)$y
    tb <- approx(x = dft, y = yb, xout = df)$y
    Tkdf[i] <- extrapolFN(tt, tb, r[i], c)
  }
  Tkdf[1] <- qnorm(0.05, lower.tail = FALSE)
  STAT <- cbind(ctr = STATISTIC)
  row.names(STAT) <- sapply(1:(k - 1), function(i) paste0("mu", 
                                                          i))
  STATCRIT <- cbind(ctr = Tkdf)
  row.names(STATCRIT) <- row.names(STAT)
  DAT <- data.frame(xold, g)
  METHOD <- c("Shirley-Williams test")
  parameter <- Inf
  names(parameter) <- "df"
  ans <- list(method = METHOD, data.name = DNAME, crit.value = STATCRIT, 
              statistic = STAT, parameter = parameter, alternative = alternative, 
              dist = "t'", model = DAT)
  class(ans) <- "osrt"
  return(ans)
}