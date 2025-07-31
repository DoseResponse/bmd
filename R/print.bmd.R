#' Print Method for bmd Objects
#'
#' Prints a summary of a bmd model object.
#'
#' @param x A bmd model object
#' @param ... Additional arguments (not used)
#' @param digits Number of significant digits to use for printing values
#'
#' @return Invisibly returns the x object
#' @export
"print.bmd" <- function(x, ..., digits = max(3, getOption("digits") - 3)) 
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
