.williamsTest <- function(x, g, alternative = c("auto", "greater", "less"), ...) 
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
  xi <- tapply(x, g, mean, na.rm = T)
  ni <- tapply(x, g, length)
  k <- nlevels(g)
  kk <- k - 1
  if (kk > 10) 
    stop("Critical t-values are only available for up to 10 dose levels.")
  N <- sum(ni)
  df <- N - k
  s2i <- tapply(x, g, var)
  s2in <- 1/df * sum(s2i * (ni - 1))
  # xiiso <- .Fortran("pava", y = as.double(xi), w = as.double(ni), 
  #                   kt = integer(k), n = as.integer(k))$y
  xiiso <- .pava(y = as.double(xi), w = as.double(ni), 
                 kt = integer(k))$y
  mui <- rep(NA, k)
  for (i in 1:k) {
    v <- k
    tmp <- rep(NA, length(1:i))
    for (u in 1:i) {
      j <- u
      tmp01 <- sapply(i:k, function(v) sum(ni[j:v] * xiiso[j:v])/sum(ni[j:v]))
      tmp[u] <- min(tmp01)
    }
    mui[i] <- max(tmp, na.rm = TRUE)
  }
  Tk <- sapply(2:k, function(i) {
    (mui[i] - xi[1])/sqrt((s2in/ni[i] + s2in/ni[1]))
  })
  extrapolFN <- function(Tki, beta, r, c) {
    out <- Tki - 0.01 * beta * (1 - r/c)
    return(out)
  }
  c <- ni[1]
  r <- ni[2:k]
  o <- c/r
  for (i in 1:kk) {
    if (o[i] < 1 | o[i] > 6) {
      warning(paste0("Ratio n0 / n", i, " is ", o[i], " and outside the range.\n\n                       Test results may not be accurate."))
    }
  }
  nrows <- nrow(williams.tk005) # PMCMRplus:::TabCrit$williams.tk005
  Tkdf <- numeric(kk)
  dft <- as.numeric(williams.tk005$rowname) # PMCMRplus:::TabCrit$williams.tk005
  xx <- c(2:6, 8, 10)
  for (i in 2:kk) {
    if (i <= 6 | i == 8 | i == 10) {
      yt <- williams.tk005[, paste0("X", i)] # PMCMRplus:::TabCrit$williams.tk005 # williams.tk005[, paste0(i)]
      yb <- williams.beta005[, paste0("X", i)] # PMCMRplus:::TabCrit$williams.beta005
    }
    else {
      yt <- sapply(1:nrows, function(j) {
        approx(x = xx, y = williams.tk005[j,-1], xout = i)$y # PMCMRplus:::TabCrit$williams.tk005
      })
      yb <- sapply(1:nrows, function(j) {
        approx(x = xx, y = williams.beta005[j,-1], xout = i)$y # PMCMRplus:::TabCrit$williams.beta005
      })
    }
    tt <- approx(x = dft, y = yt, xout = df)$y
    tb <- approx(x = dft, y = yb, xout = df)$y
    Tkdf[i] <- extrapolFN(tt, tb, r[i], c)
  }
  Tkdf[1] <- qt(0.05, df = df, lower.tail = FALSE)
  STAT <- cbind(ctr = Tk)
  row.names(STAT) <- sapply(1:(k - 1), function(i) paste0("mu", i))
  STATCRIT <- cbind(ctr = Tkdf)
  row.names(STATCRIT) <- row.names(STAT)
  parameter = c(df)
  names(parameter) <- c("df")
  METHOD <- paste("Williams trend test")
  ans <- list(method = METHOD, data.name = DNAME, crit.value = STATCRIT, 
              statistic = STAT, parameter = parameter, alternative = alternative, 
              dist = "t'")
  class(ans) <- "osrt"
  ans
}