#' Generate Bootstrap Data for BMD Analysis
#'
#' @description
#' Helper function for `bmdBoot` that generates bootstrap datasets from fitted 
#' dose-response models. Supports multiple bootstrap methods (nonparametric, 
#' parametric, and semiparametric) and handles different response types 
#' (binomial, continuous, Poisson, negative binomial).
#'
#' @param object A fitted dose-response model object (typically from `drm()`) 
#'        containing the original data and model specifications
#' @param R integer; Number of bootstrap replicates to generate (default: 1000)
#' @param bootType character; Type of bootstrap resampling method. Options are:
#'        \itemize{
#'          \item \code{"nonparametric"}: Resamples observations within dose groups
#'          \item \code{"parametric"}: Generates new data from fitted distributions
#'          \item \code{"semiparametric"}: Uses fitted values plus resampled residuals
#'        }
#' @param aggregated logical; For binomial data, whether to return aggregated 
#'        (summary) format or expanded (individual observation) format (default: TRUE)
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
#' **Count Data (Poisson/Negative Binomial):**
#' \itemize{
#'   \item Treated similarly to continuous data for nonparametric bootstrap
#'   \item Resamples within dose groups
#' }
#'
#' @section Multi-curve Support:
#' The function handles multi-curve dose-response data by preserving curve 
#' identifiers during bootstrap resampling, ensuring that the bootstrap samples 
#' maintain the original experimental structure.
#'
#' @note
#' This is an internal helper function for `bmdBoot`. It assumes the input object 
#' has the standard structure from `drm()` fitting, including components like 
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
#' # Fit a dose-response model
#' model <- drm(response/total ~ dose, weights = total, 
#'              data = binomial_data, fct = LL.4())
#' 
#' # Generate nonparametric bootstrap samples
#' boot_data <- bootDataGen(model, R = 100, bootType = "nonparametric")
#' 
#' # Generate parametric bootstrap samples  
#' boot_data_param <- bootDataGen(model, R = 100, bootType = "parametric")
#' 
#' # For continuous data with semiparametric bootstrap
#' cont_model <- drm(response ~ dose, data = continuous_data, fct = LL.4())
#' boot_data_semi <- bootDataGen(cont_model, R = 100, bootType = "semiparametric")
#' 
#' # Access first bootstrap sample
#' first_sample <- boot_data[[1]]
#' } 
bootDataGen <- function(object, R=1000, bootType="nonparametric",aggregated=TRUE){
  if(bootType=="nonparametric"){
    if(object$type=="binomial"){
      data.str <- object$data
      data.str[["number"]] <- data.str[,2]*data.str[["weights"]]
      data.e<-expandBinomial(data.str, 
                          number = "number",
                          total = "weights",
                          dose = as.character(object$call$formula[[3]]),
                          curveid = as.character(object$call$curveid))
  data.e[,"row.num"]<-1:dim(data.e)[1]
  tmp.data <- list()
  for(i in 1:R){
    sampled.expand <- data.e[as.numeric(unlist(aggregate(row.num ~ data.e[,as.character(object$call$formula[[3]])], data=data.e, 
                                                        FUN=function(x) sample(x,replace=TRUE))[[2]])),]
    if(aggregated){
      df <- aggregate(cbind(sampled.expand[,"number"],
                          sampled.expand[,"weights"]) ~ 
                      sampled.expand[,as.character(object$call$formula[[3]])],FUN = sum)
    colnames(df) <- c(as.character(object$call$formula[[3]]),
                      as.character(object$call$formula[[2]])[[2]],
                      as.character(object$call$formula[[2]])[[3]])
    tmp.data[[i]] <- df
    } else {
      if(is.null(object$call$curveid)){
        df <- data.frame(sampled.expand[,as.character(object$call$formula[[3]])],
                    sampled.expand[,"number"],
                    sampled.expand[,"weights"])
        colnames(df) <- c(as.character(object$call$formula[[3]]),
                          as.character(object$call$formula[[2]])[[2]],
                          as.character(object$call$formula[[2]])[[3]])
      } else {
        df <- data.frame(sampled.expand[,as.character(object$call$formula[[3]])],
                         sampled.expand[,"number"],
                         sampled.expand[,"weights"],
                         sampled.expand[,as.character(object$call$curveid)]
                         )
        colnames(df) <- c(as.character(object$call$formula[[3]]),
                          as.character(object$call$formula[[2]])[[2]],
                          as.character(object$call$formula[[2]])[[3]],
                          as.character(object$call$curveid))
      }
      tmp.data[[i]] <- df
    }
  }
    }
    if(object$type %in% c("continuous","Poisson","negbin1","negbin2")){
      # data.e<-object$data
      data.e<-object$origData
      data.e[,"row.num"]<-1:dim(data.e)[1]
      data.e[,"dose"]<-data.e[,as.character(object$call$formula[[3]])]
      tmp.data <- list()
      for(i in 1:R){
        tmp.data[[i]] <- data.e[as.numeric(unlist(aggregate(row.num ~ dose, data=data.e, 
                                                             FUN=function(x) sample(x,replace=TRUE))[[2]])),]
         }
    }
  } else if(bootType=="parametric"){
    if(object$type=="binomial"){
    Y <- object$data[[as.character(object$call$formula)[[2]]]]*object$data[["weights"]]
    N <- object$data[["weights"]]
    shrinks <- which(Y==N | Y==0)
    Y[shrinks] <- Y[shrinks]+0.25
    N[shrinks] <- N[shrinks]+0.5
    prob <- rep(Y/N,N)
    tmp.data <- list()
    for(i in 1:R){
      sampled.expand <- data.frame(number = rbinom(length(prob),1,prob), 
                                   dose = rep(object$data[,as.character(object$call$formula[[3]])],N), 
                                   total = 1)
      if(aggregated){
      df <- aggregateBinomial(number/total~dose, sampled.expand)
      colnames(df) <- c(as.character(object$call$formula[[3]]),
                        as.character(object$call$formula[[2]])[[2]],
                        as.character(object$call$formula[[2]])[[3]])
      tmp.data[[i]] <- df
      } else {
        colnames(sampled.expand) <- c(as.character(object$call$formula[[2]])[[2]],
                                      as.character(object$call$formula[[3]]),
                                      as.character(object$call$formula[[2]])[[3]])
        tmp.data[[i]] <- sampled.expand
      }
    }
    }
    if(object$type=="continuous"){
      origDose <- object$dataList$dose
      mean.Y <- aggregate(object$data[,as.character(object$call$formula[[2]])] ~ 
                            object$data[,as.character(object$call$formula[[3]])],
                          FUN=function(x) mean(x,na.rm=TRUE))[,2]
      sd.Y <- aggregate(object$data[,as.character(object$call$formula[[2]])] ~ 
                          object$data[,as.character(object$call$formula[[3]])],
                        FUN=function(x) sd(x,na.rm=TRUE))[,2]
      Dose<- aggregate(object$data[,as.character(object$call$formula[[2]])] ~ 
                         object$data[,as.character(object$call$formula[[3]])],
                       FUN=function(x) length(!is.na(x)))[,1]
      N.dose<- aggregate(object$data[,as.character(object$call$formula[[2]])] ~ 
                           object$data[,as.character(object$call$formula[[3]])],
                         FUN=function(x) length(!is.na(x)))[,2]
      tmp.data <- list()
      for(i in 1:R){
        sampled <- data.frame(y = rnorm(sum(N.dose),mean=rep(mean.Y,N.dose),sd=rep(sd.Y,N.dose)), 
                                     dose = rep(Dose,N.dose))
        colnames(sampled) <- c(as.character(object$call$formula[[2]]), as.character(object$call$formula[[3]]))
        tmp.data[[i]] <- sampled
        }
    }
  }
  else if(bootType=="semiparametric"){
    if(object$type=="binomial"){
      stop(paste("semiparametric is not possible for binomial data", sep=""))
    }
    if(object$type=="continuous"){
      data.st<-object$data
      
      tmp.data <- list()
      for(i in 1:R){
        sampled <- data.frame(y = fitted(object)+sample(resid(object),replace=TRUE), 
                              dose = object$data[,as.character(object$call$formula[[3]])])
        colnames(sampled) <- c(as.character(object$call$formula[[2]]), as.character(object$call$formula[[3]]))
        tmp.data[[i]] <- sampled
      }
    }
  }
tmp.data      
}


