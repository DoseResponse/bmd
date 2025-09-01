#' Generate Bootstrap Data for meta-analytic BMD Analysis
#'
#' @description
#' Helper function for `bmdBoot` that generates bootstrap datasets from fitted 
#' dose-response models. Supports multiple bootstrap methods (nonparametric, 
#' parametric, semiparametric and wild) and handles different response types 
#' (binomial, continuous).
#'
#' @param object A fitted dose-response model object of class "drcMMRE" (typically from `drmMMRE()`) 
#'        containing the original data and model specifications
#' @param R integer; Number of bootstrap replicates to generate (default: 1000)
#' @param bootType character; Type of bootstrap resampling method. Options are:
#'        \itemize{
#'          \item \code{"nonparametric"}: Resamples observations within dose groups
#'          \item \code{"parametric"}: Generates new data from fitted distributions
#'          \item \code{"semiparametric"}: Uses fitted values plus resampled residuals
#'          \item \code{"wild"}: Uses fitted values plus resampled residuals
#'        }
#'
#' @details
#' The function implements different bootstrap strategies based on data type and method:
#' 
#' **Nonparametric Bootstrap:**
#' \itemize{
#'   \item \strong{Binomial}: Expands binomial data to individual observations, 
#'         resamples within dose groups, then optionally re-aggregates
#'   \item \strong{Continuous/Count}: Resamples observations within dose groups 
#'         from the original dataset
#' }
#' 
#' **Parametric Bootstrap:**
#' \itemize{
#'   \item \strong{Binomial}: Generates new binomial observations using estimated 
#'         success probabilities, with continuity correction (adds 0.25/0.5) for 
#'         boundary cases
#'   \item \strong{Continuous}: Generates normal random variables using 
#'         dose-specific means and standard deviations from original data
#' }
#' 
#' **Semiparametric Bootstrap:**
#' \itemize{
#'   \item \strong{Continuous only}: Uses fitted values plus resampled residuals
#'   \item \strong{Binomial}: Not supported (throws error)
#' }
#'
#' @return A list of length R containing bootstrap datasets. Each element is a 
#'         data.frame with the same structure as the original data, containing:
#'         \itemize{
#'           \item Dose variable (same name as original)
#'           \item Response variable(s) (same name(s) as original)
#'           \item For binomial data: number of successes and total observations
#'           \item For multi-curve data: curve identifier (if present)
#'         }
#'
#' @section Data Type Handling:
#' 
#' **Binomial Data:**
#' \itemize{
#'   \item Handles both aggregated (n successes out of N trials) and expanded formats
#'   \item Preserves dose group structure during resampling
#'   \item Applies continuity correction in parametric bootstrap
#' }
#' 
#' **Continuous Data:**
#' \itemize{
#'   \item Maintains dose group structure
#'   \item Preserves within-group variability patterns
#'   \item Uses original data (`origData`) when available
#' }
#' 
#'
#' @note
#' This is an internal helper function for `bmdBoot`. It assumes the input object 
#' has the standard structure from `drmMMRE()` fitting, including components like 
#' `data`, `origData`, `call`, `type`, etc.
#' 
#' **Important considerations:**
#' \itemize{
#'   \item Semiparametric bootstrap requires model residuals and fitted values
#'   \item Parametric bootstrap assumes distributional assumptions are met
#'   \item Large R values may require substantial memory for complex datasets
#' }
#'
#' @seealso 
#' \code{\link{bmdBoot}} for the main bootstrap BMD function,
#'
#' @examples
#' \dontrun{
#' # Typically called internally by bmdBoot, but can be used directly:
#' 
#' set.seed(1)
#' data0 <- data.frame(x = rep(drcData::ryegrass$conc, 2),
#'                     y = rep(drcData::ryegrass$rootl, 2) +
#'                       c(rnorm(n = nrow(drcData::ryegrass), mean = 2, sd = 0.5),
#'                         rnorm(n = nrow(drcData::ryegrass), mean = 2.7, sd = 0.7)),
#'                     EXP_ID = rep(as.character(1:2), each = nrow(drcData::ryegrass)))
#' 
#' modMMRE <- drmMMRE(y~x, exp_id = EXP_ID, data = data0, fct = LL.4())
#' boot_data <- bootDataGenMMRE(modMMRE, R = 1000, bootType = "parametric")
#' 
#' # Access first bootstrap sample
#' first_sample <- boot_data[[1]]
#' } 
bootDataGenMMRE <- function(object, R=1000, bootType=c("nonparametric", "semiparametric", "parametric", "wild")){
  if(!inherits(object, "drcMMRE")){
    stop("bootDataGenMMRE only works for object of type drcMMRE")
  }
  
  R <- as.integer(R)
  if(is.na(R)){
    stop("R must be a positive integer")
  }
  if(R<=0){
    stop("R must be a positive integer")
  }
  
  bootType <- match.arg(bootType)
  
  # resample for each experiment independently, then collect individual resampled datasets
  boot.data.per.exp_id <- lapply(names(object$objList), 
                                 function(x){
                                   object0 <- object$objList[[x]]
                                   object0$call$formula <- object$call$formula
                                   boot.data.per.exp_id0 <- bootDataGen(object0,R=R, bootType=bootType, aggregated = FALSE)
                                   # add exp_id to data
                                   boot.data.per.exp_id0 <- lapply(boot.data.per.exp_id0,
                                                                   function(z) {
                                                                     z[[as.character(object$call$exp_id)]] <- rep(x, nrow(z))
                                                                     z
                                                                     }
                                   )
                                   boot.data.per.exp_id0
                                   }
                                 )
  # collect individual resampled datasets
  tmp.data <- lapply(
    1:R,
    function(r) do.call(rbind, lapply(boot.data.per.exp_id, function(list) list[[r]]))
    )
  
  tmp.data
}


