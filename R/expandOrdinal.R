expandOrdinal <- function(object){
  df.small <- object$data
  df.small$row.orig <- 1:nrow(df.small)
  df.long <- reshape2:::melt(df.small, measure.vars = object$categories)
  idx <- rep(1:nrow(df.long), df.long$value)
  df.long.expand <- df.long[idx,] 
  df.long.expand
}
