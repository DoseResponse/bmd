#' Extract Coefficients from Heteroscedastic Dose-Response Models
#'
#' @description
#' S3 method to extract coefficients from fitted heteroscedastic dose-response models.
#' This method extends the generic `coef()` function to handle `drcHetVar` objects,
#' which contain both dose-response curve parameters and variance model parameters.
#'
#' @param object An object of class "drcHetVar" containing a fitted heteroscedastic 
#'        dose-response model with separate curve and variance parameters
#' @param ... Additional arguments (currently unused but included for S3 method consistency)
#'
#' @details
#' Heteroscedastic dose-response models fit two components simultaneously:
#' \itemize{
#'   \item \strong{Curve parameters}: Describe the dose-response relationship 
#'         (e.g., slope, inflection point, asymptotes)
#'   \item \strong{Variance parameters}: Model how the response variance changes 
#'         with dose or fitted values
#' }
#' 
#' This method combines both parameter sets into a single named vector for easy 
#' access and interpretation. Parameter names are prefixed to distinguish between 
#' the two model components:
#' \itemize{
#'   \item \code{"Curve:"} prefix for dose-response curve parameters
#'   \item \code{"Sigma:"} prefix for variance model parameters
#' }
#'
#' @return A named numeric vector containing all model parameters:
#' \itemize{
#'   \item Curve parameters with names prefixed by "Curve:"
#'   \item Variance parameters with names prefixed by "Sigma:"
#' }
#' 
#' @note
#' This method assumes the `drcHetVar` object contains:
#' \itemize{
#'   \item \code{curvePar}: Named vector of curve parameters
#'   \item \code{sigmaPar}: Named vector of variance parameters
#' }
#' 
#' The original parameter names from the fitting procedure are preserved but 
#' prefixed for clarity.
#'
#' @seealso 
#' \code{\link{coef}} for the generic coefficient extraction function,
#'
#' @examples
#' \dontrun{
#' # Assuming you have a fitted heteroscedastic model
#' }
#'
#' @method coef drcHetVar
#' @export
coef.drcHetVar <- function(object, ...){
  curvePar <- object$curvePar
  names(curvePar) <- paste0("Curve:", names(curvePar))
  sigmaPar <- object$sigmaPar
  names(sigmaPar) <- paste0("Sigma:", names(sigmaPar))
  c(curvePar, sigmaPar)
}
