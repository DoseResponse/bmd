.pava <- function(y, w, kt) {
  n <- length(y)
  
  # Initialize kt
  kt <- 1:n
  
  if (n > 1) {
    for (i in 2:n) {
      if (y[i - 1] > y[i]) {
        k1 <- kt[i]
        k2 <- kt[i - 1]
        
        # Update kt
        for (j in 1:n) {
          if (kt[j] == k1) {
            kt[j] <- k2
          }
        }
        
        # Update y and w
        wnew <- w[i - 1] + w[i]
        ynew <- (w[i - 1] * y[i - 1] + w[i] * y[i]) / wnew
        for (j in 1:n) {
          if (kt[j] == k2) {
            y[j] <- ynew
            w[j] <- wnew
          }
        }
      }
    }
  }
  
  return(list(y = y, w = w, kt = kt))
}
