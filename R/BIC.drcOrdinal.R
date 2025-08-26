#' BIC Method for drcOrdinal Objects
#'
#' @description
#' S3 method to calculate the Bayesian Information Criterion (BIC) for fitted 
#' ordinal dose-response models. This method extends the generic `BIC()` function 
#' to handle `drcOrdinal` objects which contain multiple fitted dose-response models 
#' for different response categories.
#'
#' @param object An object of class "drcOrdinal" containing fitted ordinal 
#'        dose-response models
#' @param ... Additional arguments passed to the method. Currently supports:
#'        \itemize{
#'          \item \code{epsilon}: Small positive value to avoid log(0) in likelihood 
#'                calculations (default: 1e-16)
#'        }
#'
#' @details
#' The BIC is calculated using the standard formula:
#' \deqn{BIC = k \log(n) - 2 \log(L)}{BIC = k * log(n) - 2 * log(L)}
#' where:
#' \itemize{
#'   \item \eqn{k}{k} is the total number of parameters across all models in the ordinal fit
#'   \item \eqn{n}{n} is the total number of observations (sum of weights)
#'   \item \eqn{L}{L} is the likelihood of the fitted model
#' }
#' 
#' For ordinal dose-response models, the total number of parameters is the sum 
#' of parameters from all individual models in the `drmList` component, and the 
#' sample size is calculated as the sum of weights in the data.
#' 
#' The `epsilon` parameter is used in the `logLik()` calculation to prevent 
#' numerical issues when probabilities approach zero.
#'
#' @return A numeric value representing the BIC of the fitted ordinal dose-response model.
#'         Lower values indicate better model fit when comparing models.
#'
#' @section Model Comparison:
#' BIC penalizes model complexity more heavily than AIC, making it useful for:
#' \itemize{
#'   \item Comparing ordinal dose-response models with different numbers of parameters
#'   \item Model selection in situations where simpler models are preferred
#'   \item Avoiding overfitting in dose-response modeling
#' }
#'
#' @note
#' This method assumes that the `drcOrdinal` object contains:
#' \itemize{
#'   \item A `drmList` component with fitted dose-response models
#'   \item A `data` component with the original data including weights
#'   \item A `weights` component specifying the column name for observation weights
#' }
#'
#' @seealso 
#' \code{\link{BIC}} for the generic BIC function,
#' \code{\link{AIC.drcOrdinal}} for the corresponding AIC method,
#' \code{\link{logLik.drcOrdinal}} for the log-likelihood method
#'
#' @examples
#' \dontrun{
#' # Assuming you have a fitted drcOrdinal object
#' ordinal_model <- drm(response ~ dose, 
#'                      weights = count,
#'                      data = ordinal_data, 
#'                      fct = cumulative.logit())
#' 
#' # Calculate BIC
#' bic_value <- BIC(ordinal_model)
#' 
#' # Compare models
#' model1_bic <- BIC(ordinal_model1)
#' model2_bic <- BIC(ordinal_model2)
#' 
#' # Lower BIC indicates better model
#' if(model1_bic < model2_bic) {
#'   cat("Model 1 is preferred\n")
#' } else {
#'   cat("Model 2 is preferred\n")
#' }
#' 
#' # Use custom epsilon for numerical stability
#' bic_custom <- BIC(ordinal_model, epsilon = 1e-12)
#' }
#'
#' @method BIC drcOrdinal
#' @export
BIC.drcOrdinal <- function(object, ...){
  dots <- list(...)
  if (!is.null(dots$epsilon)){
    epsilon <- dots$epsilon
  } else {
    epsilon <- 1e-16
  }
  
  n.parameters <- sum(sapply(object$drmList, function(mod) length(mod$coefficients)))
  n.obs <- sum(object$data[[object$weights]])
  n.parameters * log(n.obs) - 2 * logLik(object, epsilon) 
}
