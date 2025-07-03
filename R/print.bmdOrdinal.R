
#' @title S3 method
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
