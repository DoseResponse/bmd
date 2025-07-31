#' Generate Bootstrap Data for Ordinal Dose-Response Analysis
#'
#' @description
#' Helper function for `bmdOrdinal` and `bmdOrdinalMA` that generates bootstrap 
#' datasets from fitted ordinal dose-response models. Supports multiple bootstrap 
#' methods tailored for ordinal response data with different resampling strategies.
#'
#' @param object A fitted ordinal dose-response model object containing the original 
#'        data, model specifications, response levels, and fitted parameters
#' @param R integer; Number of bootstrap replicates to generate (default: 500)
#' @param bootType character; Type of bootstrap resampling method. Options are:
#'        \itemize{
#'          \item \code{"nonparametric"}: Resamples individual observations within dose groups
#'          \item \code{"parametric"}: Generates new observations using estimated category probabilities
#'          \item \code{"model"}: Generates observations using fitted model probabilities
#'          \item \code{"hierarchical"}: Hierarchical resampling accounting for block structure
#'        }
#'
#' @details
#' The function implements four distinct bootstrap strategies for ordinal data:
#' 
#' **Nonparametric Bootstrap:**
#' \itemize{
#'   \item Expands ordinal data to individual observations
#'   \item Resamples observations within each original data row
#'   \item Maintains the empirical distribution within dose groups
#'   \item Most conservative approach, makes no distributional assumptions
#' }
#' 
#' **Parametric Bootstrap:**
#' \itemize{
#'   \item Estimates category probabilities from original data
#'   \item Applies continuity correction for boundary cases:
#'     \deqn{p_{corrected} = \frac{1/K^2}{n + 1/K} \text{ when } p = 0}{p_corrected = (1/K^2)/(n + 1/K) when p = 0}
#'     \deqn{p_{corrected} = \frac{n + 1/K^2}{n + 1/K} \text{ when } p = 1}{p_corrected = (n + 1/K^2)/(n + 1/K) when p = 1}
#'   where K is the number of response categories and n is the sample size
#'   \item Generates new observations using multinomial sampling
#' }
#' 
#' **Model-Based Bootstrap:**
#' \itemize{
#'   \item Uses fitted model probabilities (\code{object$pFun}) for each dose
#'   \item Generates observations directly from the fitted dose-response relationship
#'   \item Assumes the fitted model accurately represents the true relationship
#' }
#' 
#' **Hierarchical Bootstrap:**
#' \itemize{
#'   \item Accounts for block/cluster structure in the data
#'   \item Requires \code{object$blocks} to be specified
#'   \item Uses weighted resampling within dose groups
#'   \item Maintains hierarchical correlation structure
#' }
#'
#' @return A list of length R containing bootstrap datasets. Each element is a 
#'         data.frame with the same structure as the original ordinal data, containing:
#'         \itemize{
#'           \item Dose variable (same name as original)
#'           \item Count columns for each response category
#'           \item For hierarchical: block identifier and total counts
#'           \item All columns from original data structure preserved
#'         }
#'
#' @section Data Processing:
#' The function uses several data manipulation steps:
#' \itemize{
#'   \item \strong{Expansion}: Converts aggregated ordinal data to individual observations
#'   \item \strong{Resampling}: Applies the specified bootstrap method
#'   \item \strong{Aggregation}: Converts back to count format using \code{reshape2::dcast}
#'   \item \strong{Completion}: Ensures all response categories are present (fills with 0 if missing)
#' }
#'
#' @section Dependencies:
#' This function requires the following packages:
#' \itemize{
#'   \item \code{reshape2}: For data reshaping operations
#'   \item \code{dplyr}: For data manipulation (hierarchical bootstrap only)
#'   \item \code{tidyr}: For data tidying operations (hierarchical bootstrap only)
#' }
#' 
#' The function will stop with an informative error if required packages are not installed.
#'
#' @section Continuity Correction:
#' For parametric bootstrap, when category probabilities are 0 or 1, a continuity 
#' correction is applied to prevent degenerate sampling. This ensures all categories 
#' have some probability of being selected, improving bootstrap stability.
#'
#' @note
#' This is an internal helper function for ordinal BMD bootstrap procedures. It assumes 
#' the input object has the standard structure from ordinal dose-response fitting, 
#' including components like \code{data}, \code{levels}, \code{dose}, \code{pFun}, etc.
#' 
#' **Method Selection Guidelines:**
#' \itemize{
#'   \item Use \code{"nonparametric"} for robust, assumption-free bootstrap
#'   \item Use \code{"parametric"} when sample sizes are small
#'   \item Use \code{"model"} to assess model-based uncertainty
#'   \item Use \code{"hierarchical"} for clustered/blocked experimental designs
#' }
#'
#' @seealso 
#' \code{\link{bmdOrdinal}} for ordinal BMD estimation,
#' \code{\link{bmdOrdinalMA}} for model-averaged ordinal BMD,,
#' \code{\link[reshape2]{dcast}} for data reshaping
#'
#' @examples
#' \dontrun{
#' # Typically called internally, but can be used directly:
#' 
#' # Assume you have a fitted ordinal model
#' ordinal_model <- fitOrdinalModel(response ~ dose, data = ordinal_data)
#' 
#' # Generate nonparametric bootstrap samples
#' boot_data_np <- bootDataGenOrdinal(ordinal_model, R = 100, 
#'                                    bootType = "nonparametric")
#' 
#' # Generate parametric bootstrap samples
#' boot_data_param <- bootDataGenOrdinal(ordinal_model, R = 100, 
#'                                       bootType = "parametric")
#' 
#' # Generate model-based bootstrap samples
#' boot_data_model <- bootDataGenOrdinal(ordinal_model, R = 100, 
#'                                       bootType = "model")
#' 
#' # For hierarchical data with blocks
#' boot_data_hier <- bootDataGenOrdinal(ordinal_model, R = 100, 
#'                                      bootType = "hierarchical")
#' 
#' # Access first bootstrap sample
#' first_sample <- boot_data_np[[1]]
#' head(first_sample)
#' }
bootDataGenOrdinal <- function(object, R = 500, bootType = c("nonparametric", "parametric", "model", "hierarchical")){
  bootType <- match.arg(bootType)
  ## avoiding global binding of variables issues
  variable <- NULL
  row.num <- NULL
  value <- NULL
  row.orig <- NULL
  rm(list=c("variable", "row.num", "value", "row.orig"))
  if(!requireNamespace("reshape2")){
    stop('package "reshape2" must be installed to use bootstrapping with ordinal dose-response model')
  }
  if(!requireNamespace("dplyr")){
    stop('package "dplyr" must be installed to use bootstrapping with ordinal dose-response model')
  }
  if(!requireNamespace("tidyr")){
    stop('package "tidyr" must be installed to use bootstrapping with ordinal dose-response model')
  }
  
  if (bootType == "nonparametric") {
    data.e <- expandOrdinal(object)
    data.e[, "row.num"] <- 1:nrow(data.e)
    tmp.data <- list()
    for (i in 1:R) {
      sampled.expand <- data.e[as.numeric(unlist(aggregate(row.num ~ 
                                                             data.e[, "row.orig"], 
                                                           data = data.e, 
                                                           FUN = function(x) sample(x, 
                                                                                    replace = TRUE))[[2]])), ]
      columns.rem <- colnames(subset(sampled.expand, select=-c(variable,row.num,value,row.orig)))
      df <- reshape2::dcast(sampled.expand, as.formula(paste(paste(columns.rem, collapse = "+")," ~ variable")), length)
      for(j in 1:length(object$levels)){
        if(!(object$levels[j] %in% colnames(df))){
          df[,object$levels[j]]<-0
        }
      }
      tmp.data[[i]] <- df
    }
  }
  if (bootType == "parametric") {
    data.e <- expandOrdinal(object)
    #data.e[, "row.num"] <- 1:dim(data.e)[1]
    tmp.data <- list()
    for (i in 1:R) {
      p0 <- aggregate(variable ~ data.e[, "row.orig"],
                      data = data.e, FUN = function(x) table(x)/length(x))
      sampled.expand <- data.e
      for(j in 1:length(unique(data.e$row.orig))){
        data.size <- length(sampled.expand$variable[data.e$row.orig==j])
        prop0 <- p0[j,-1]
        prop0[prop0==0] <- (1/length(object$levels)^2)/(data.size+1/length(object$levels)) # (1/4)/(data.size + 1/2)
        prop0[prop0==1] <- (data.size + 1/length(object$levels)^2)/(data.size+1/length(object$levels)) # (data.size + 1/4)/(data.size+1/2)
        sampled.expand$variable[sampled.expand$row.orig==j] <-
          #rep(unlist(object$levels), as.numeric(rmultinom(1, size = data.size, prob = prop0)))
          unlist(sample(object$levels, size = data.size, replace = TRUE, prob = prop0))
      }
      columns.rem <- colnames(subset(sampled.expand, select=-c(variable,value,row.orig)))
      df <- reshape2::dcast(sampled.expand, as.formula(paste(paste(columns.rem, collapse = "+")," ~ variable")), length)
      for(j in 1:length(object$levels)){
        if(!(object$levels[j] %in% colnames(df))){
          df[,object$levels[j]]<-0
        }
      }
      tmp.data[[i]] <- df
    }
  }
  if (bootType == "model") {
    data.e <- expandOrdinal(object)
    tmp.data <- list()
    for (i in 1:R) {
      sampled.expand <- data.e
      sampled.expand[, "variable"] <- sapply(data.e[,object$dose], function(x) unlist(sample(object$levels, size = 1, replace = TRUE, object$pFun(x))))
      columns.rem <- colnames(subset(sampled.expand, select=-c(variable,value,row.orig)))
      df <- reshape2::dcast(sampled.expand, as.formula(paste(paste(columns.rem, collapse = "+")," ~ variable")), length)
      for(j in 1:length(object$levels)){
        if(!(object$levels[j] %in% colnames(df))){
          df[,object$levels[j]]<-0
        }
      }
      tmp.data[[i]] <- df
    }
  }
  if(bootType == "hierarchical"){
    if(is.null(object$blocks)){
      stop('ordinal dose-response model does not include blocks. Hierarchical resampling is not possible.')
      # cat("\n", 'Argument "block" needs to be specified for hierarchical resampling.', "\n")
      # tmp.data <- NULL
    } 
    
    resample_fun <- function(levels, dose, weights, blocks, data){
      name <- NULL
      rm(list("name"))
      data %>%
        dplyr::mutate(row.orig = 1:n()) %>% 
        dplyr::group_by(.data[[dose]]) %>% 
        dplyr::slice_sample(prop = 1, weight_by = eval(parse(text=weights)), replace = TRUE) %>% 
        tidyr::pivot_longer(levels) %>% 
        dplyr::group_by(row.orig) %>% 
        tidyr::uncount(value) %>% 
        dplyr::group_by(row.orig) %>% 
        dplyr::slice_sample(prop = 1, replace = TRUE) %>% 
        dplyr::group_by(.data[[dose]], .data[[weights]], .data[[blocks]], row.orig) %>% 
        dplyr::count(name) %>% 
        dplyr::mutate(total = sum(n)) %>% 
        tidyr::pivot_wider(names_from = name, values_from = n, values_fill = 0)
    }
    
    tmp.data <- list()
    for (i in 1:R){
      tmp.data[[i]] <- resample_fun(object$levels, object$dose, object$weights, object$blocks, object$data)
    }
  }
  tmp.data
}
