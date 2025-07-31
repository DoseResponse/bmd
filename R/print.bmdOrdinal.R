#' Print Method for bmdOrdinalr Objects
#'
#' Prints a summary of a bmdOrdinal model object.
#'
#' @param x A bmdOrdinal model object
#' @param ... Additional arguments (not used)
#' @param digits Number of significant digits to use for printing values
#'
#' @return Invisibly returns the x object
#' @export
print.bmdOrdinal <- function(x, ..., digits = max(3, getOption("digits") - 3)) 
{
  object <- x
  classList <- class(object)
  
  if (length(object$Results)>0) 
  {
    cat("\n")
    print(object$Results)
  } else {
    cat("Problem occured. Please check whether the choice of bmr is meaningful\n")
  }
  cat("\n")
  
  invisible(object)
}
