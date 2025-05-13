print.drcMMRE <- function (x, ..., digits = max(3, getOption("digits") - 3)) 
{
  object <- x
  classList <- class(object)
  cat(paste("\n", "A 'drcMMRE' model.", "\n", sep = ""))
  cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
  if (length(coef(object)) > 0) {
    cat("Coefficients:\n")
    print.default(format(coef(object), digits = digits), 
                  print.gap = 2, quote = FALSE)
  }
  else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(object)
}