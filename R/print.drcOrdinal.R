print.drcOrdinal <- function(object, ..., digits = max(3, getOption("digits") - 3)) 
{
  classList <- class(object)
  cat("\n", "A 'drcOrdinal' model.", "\n", 
      "Categories: ", paste(object$categories, collapse = ", "), "\n",
      "Dose: ", object$dose, "\n", 
      "Weights: ", object$weights, "\n",
      "Function: ", object$drmList[[1]]$fct$name, "\n\n",
      "This 'drcOrdinal' model is composed of the following 'drc' models:", sep="")
  
  lapply(1:length(object$drmList),
         function(i){
           cat("\n", object$catMerged[[i]], ":\n", "Coefficients:\n", sep = "")
           print.default(format(coef(object$drmList[[i]]), digits = digits), print.gap = 2, quote = FALSE)
         })
  
  invisible(object)
}
