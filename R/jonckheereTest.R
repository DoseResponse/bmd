.jonckheereTest <- function (x, g, alternative = c("auto", "two.sided", "greater", "less"), 
          continuity = FALSE, ...) 
{
  # Prepare observations
  if (length(x) != length(g)) 
    stop("'x' and 'g' must have the same length")
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  OK <- complete.cases(x, g)
  x <- x[OK]
  g <- g[OK]
  if (!all(is.finite(g))) 
    stop("all group levels must be finite")
  
  # Find trend direction
  alternative <- match.arg(alternative)
  if(alternative == "auto"){
    lm0 <- lm(x ~ g)$coef[["g"]]
    alternative <- ifelse(lm0 > 0, "greater", "less")
  }
  
  g <- factor(g)
  k <- nlevels(g)
  if (k < 2) 
    stop("all observations are in the same group")
  if (!is.logical(continuity)) 
    stop("'continuity' must be 'FALSE' or 'TRUE'")
  # alternative <- match.arg(alternative)
  n <- length(x)
  if (n < 2) 
    stop("needs at least 3 observations")
  o <- order(g)
  g <- g[o]
  x <- x[o]
  nij <- tapply(x, g, length)
  X <- matrix(NA, ncol = k, nrow = max(nij))
  j <- 0
  for (i in 1:k) {
    for (l in 1:nij[i]) {
      j = j + 1
      X[l, i] <- x[j]
    }
  }
  psi.f <- function(u) {
    psi <- (sign(u) + 1)/2
    psi
  }
  Uij <- function(i, j, X) {
    ni <- nij[i]
    nj <- nij[j]
    sumUij <- 0
    for (s in (1:ni)) {
      for (t in (1:nj)) {
        sumUij <- sumUij + psi.f(X[t, j] - X[s, i])
      }
    }
    sumUij
  }
  J <- 0
  for (i in (1:(k - 1))) {
    for (j in ((i + 1):k)) {
      J = J + Uij(i, j, X)
    }
  }
  mu <- (n^2 - sum(nij^2))/4
  st <- 0
  for (i in (1:k)) {
    st <- st + nij[i]^2 * (2 * nij[i] + 3)
    st
  }
  TIES <- FALSE
  TIES <- (sum(table(rank(x)) - 1) > 0)
  if (!TIES) {
    s <- sqrt((n^2 * (2 * n + 3) - st)/72)
    S <- J - mu
  } else {
    # warning("Ties are present. Jonckheere z was corrected for ties.")
    S <- J - mu
    nt <- as.vector(table(x))
    s <- sqrt((n * (n - 1) * (2 * n + 5) - sum(nij * (nij - 
                                                        1) * (2 * nij + 5)) - sum(nt * (nt - 1) * (2 * nt + 
                                                                                                     5)))/72 + (sum(nij * (nij - 1) * (nij - 2)) * sum(nt * 
                                                                                                                                                         (nt - 1) * (nt - 2)))/(36 * n * (n - 1) * (n - 2)) + 
                (sum(nij * (nij - 1)) * sum(nt * (nt - 1)))/(8 * 
                                                               n * (n - 1)))
  }
  if (continuity) {
    S <- sign(S) * (abs(S) - 0.5)
  }
  STATISTIC <- S/s
  if (alternative == "two.sided") {
    PVAL <- 2 * min(pnorm(abs(STATISTIC), lower.tail = FALSE), 
                    0.5)
  } else if (alternative == "greater") {
    PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
  } else {
    PVAL <- pnorm(STATISTIC)
  }
  ESTIMATES <- J
  names(ESTIMATES) <- "JT"
  names(STATISTIC) <- "z"
  RVAL <- list(statistic = STATISTIC, p.value = PVAL, method = "Jonckheere-Terpstra test", 
               data.name = DNAME, alternative = alternative, estimates = ESTIMATES)
  return(RVAL)
}
