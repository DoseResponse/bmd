print.drcHetVar <- function(x, ..., digits = max(3, getOption("digits") - 3)){
  object <- x
  classList <- class(object)
  cat(paste("\n", "A 'drcHetVar' model.", "\n", sep = ""))
  cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
  
  # Coeffecients
  cat("Curve Coefficients:\n")
  print.default(format(object$curvePar, digits = digits), 
                print.gap = 2, quote = FALSE)
  cat("\n")
  
  cat("Variance Coefficients:\n")
  print.default(format(object$sigmaPar, digits = digits), 
                print.gap = 2, quote = FALSE)
  cat("\n")
  
  # END
  invisible(object)
}