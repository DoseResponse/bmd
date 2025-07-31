#' Inverse Regression for Benchmark Dose Calculation
#'
#' @description
#' Helper function for `bmd` that calculates Benchmark Dose (BMD) and confidence 
#' intervals using inverse regression. This function handles multiple BMD definitions 
#' and background response types by dynamically modifying the dose-response function 
#' and its derivatives to solve for the dose corresponding to a specified benchmark response.
#'
#' @param object A fitted dose-response model object containing the model function,
#'        derivatives, coefficients, and data information
#' @param bmr numeric; Benchmark response level for which to calculate the BMD
#' @param level numeric; Confidence level for BMD confidence intervals (default: 0.9)
#' @param slope character; Direction of dose-response relationship, either 
#'        "increasing" or "decreasing"
#' @param backgType character; Type of background response calculation. Options include:
#'        "modelBased", "hybridSD", "absolute"
#' @param backg numeric; Background response value (used with certain backgType options).
#'        Default: NA
#' @param catLev numeric; Category level for multi-level responses. Default: NA
#' @param extFactor numeric; Extension factor for maximum dose search range (default: 10)
#' @param def character; BMD definition type. Options include:
#'        \itemize{
#'          \item "additional"/"added": Additional risk over background
#'          \item "excess": Excess risk over background  
#'          \item "relative": Relative change from background
#'          \item "extra": Extra risk (between background and maximum)
#'          \item "hybridAdd": Hybrid additional risk using statistical criteria
#'          \item "hybridExc": Hybrid excess risk using statistical criteria
#'          \item "point": Point estimate (not supported for inverse regression)
#'        }
#' @param useSD numeric; Standard deviation value used in hybrid definitions
#' @param sandwich.vcov logical; Whether to use sandwich variance-covariance matrix
#'        for robust standard errors (default: FALSE)
#'
#' @details
#' The function implements inverse regression by solving:
#' \deqn{f(BMD) = BMR}{f(BMD) = BMR}
#' where f is the transformed dose-response function based on the specified definition.
#' 
#' **BMD Definitions:**
#' 
#' **Additional Risk:** \deqn{f(x) = |f_0(0) - f_0(x)|}{f(x) = |f_0(0) - f_0(x)|}
#' 
#' **Excess Risk:** \deqn{f(x) = \frac{|f_0(0) - f_0(x)|}{f_0(0)}}{f(x) = |f_0(0) - f_0(x)|/f_0(0)} (decreasing)
#' \deqn{f(x) = \frac{f_0(x) - f_0(0)}{1 - f_0(0)}}{f(x) = (f_0(x) - f_0(0))/(1 - f_0(0))} (increasing)
#' 
#' **Relative Risk:** \deqn{f(x) = \frac{|f_0(0) - f_0(x)|}{f_0(0)}}{f(x) = |f_0(0) - f_0(x)|/f_0(0)}
#' 
#' **Extra Risk:** \deqn{f(x) = \frac{|f_0(0) - f_0(x)|}{|f_0(0) - f_0(\infty)|}}{f(x) = |f_0(0) - f_0(x)|/|f_0(0) - f_0(infinity)|}
#' 
#' **Hybrid Definitions:** Use normal distribution probabilities:
#' \deqn{f(x) = \Phi\left(\frac{f_0(0) - k \cdot SD - f_0(x)}{SD}\right) - (1-\Phi(k))}{f(x) = Phi((f_0(0) - k*SD - f_0(x))/SD) - (1-Phi(k))}
#' where k is the number of standard deviations (default: 2).
#'
#' @return A matrix with 1 row and 3 columns:
#' \itemize{
#'   \item \strong{BMD}: Benchmark dose estimate
#'   \item \strong{BMDL}: Lower confidence limit for BMD
#'   \item \strong{BMDU}: Upper confidence limit for BMD
#' }
#' 
#' If root finding fails for confidence limits, BMDL is set to 0 and BMDU to Inf.
#'
#' @section Function Modification:
#' The function dynamically modifies the model's derivative function to account 
#' for different BMD definitions by:
#' \itemize{
#'   \item Parsing the derivative function body
#'   \item Adding appropriate transformation terms
#'   \item Handling parameter indexing for fixed vs. estimated parameters
#'   \item Supporting both 4-parameter and 5-parameter models
#' }
#'
#' @section Confidence Intervals:
#' Confidence intervals are calculated using the delta method:
#' \deqn{CI = BMD \pm t_{\alpha/2} \sqrt{h(BMD)^T \Sigma h(BMD)}}{CI = BMD +/- t * sqrt(h(BMD)' * Sigma * h(BMD))}
#' where h(BMD) is the gradient vector and Σ is the variance-covariance matrix.
#' 
#' For continuous data, t-distribution quantiles are used; for other data types, 
#' normal distribution quantiles are used.
#'
#' @section Error Handling:
#' The function includes several error checks:
#' \itemize{
#'   \item Stops if "point" definition is used (not supported)
#'   \item Stops if hybrid definitions are used without appropriate backgType
#'   \item Stops if no background value is provided when required
#'   \item Stops if no solution exists for the specified BMR
#' }
#'
#' @note
#' This is an internal helper function for BMD calculation. It assumes the input 
#' object has the standard structure from dose-response model fitting, including 
#' components like `fct`, `dataList`, `type`, etc.
#' 
#' **Important considerations:**
#' \itemize{
#'   \item The search range is limited to `[0, extFactor × max(dose)]`
#'   \item Hybrid definitions require careful specification of useSD and backg
#'   \item Function modification uses string manipulation and may be sensitive to model structure
#' }
#'
#' @seealso 
#' \code{\link{bmd}} for the main BMD calculation function,
#' \code{\link{uniroot}} for the root finding algorithm,
#' \code{\link{vcov}} for variance-covariance matrix extraction,
#' \code{\link[sandwich]{sandwich}} for robust variance estimation
#'
#' @examples
#' \dontrun{
#' # Typically called internally by bmd(), but can be used directly:
#' # Need to be updated
#' # Fit a dose-response model
#' model <- drm(response ~ dose, data = dose_data, fct = LL.4())
#' 
#' # Calculate BMD for 10% additional risk
#' bmd_result <- invBmd(object = model, 
#'                      bmr = 0.1, 
#'                      level = 0.9,
#'                      slope = "increasing",
#'                      def = "additional",
#'                      useSD = sd(model$data$response))
#' 
#' # Extract results
#' bmd_estimate <- bmd_result[1, "BMD"]
#' bmd_lower <- bmd_result[1, "BMDL"] 
#' bmd_upper <- bmd_result[1, "BMDU"]
#' 
#' # Using hybrid definition with 1.5 standard deviations
#' bmd_hybrid <- invBmd(object = model,
#'                      bmr = 0.05,
#'                      def = "hybridAdd",
#'                      backgType = "hybridSD",
#'                      backg = 1.5,
#'                      useSD = sd(model$data$response))
#' }
invBmd <- function(object, bmr, level=0.9, slope, backgType="modelBased", 
                     backg=NA, catLev=NA, extFactor=10, def, useSD=useSD, sandwich.vcov=FALSE){
  
  ParmVec0 <- object$fct$fixed
  ParmVec <- ParmVec0
  ParmVec[is.na(ParmVec0)] <- coef(object)
  g<-object$fct$fct
  h<-object$fct$deriv1
  if(substr(deparse(as.list(body(h))[[length(as.list(body(h)))]])[1],1,6)=="return"){
    element.num <- length(as.list(body(h)))-1
    dpgh <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
    pgh <- paste(dpgh, "[, notFixed]", sep="")
    body(h)[[element.num]] <- as.call(str2lang(pgh))
    body(h)[[element.num+1]] <- substitute(derMat)
  } else {
    element.num <- length(as.list(body(h)))
  }
  
  
  if(def %in% c("additional", "added")){
    if(identical(slope,"decreasing")){
      g0 <- function(x){as.numeric(g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                     g(x, matrix(coef(object), 1, length(coef(object)))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " + as.matrix(cbind(0, 0, -1, 0))[, notFixed]", sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " + as.matrix(cbind(0, 0, -1, 0, 0))[, notFixed]", sep="")
      }
      body(h)[[element.num]] <- as.call(str2lang(pgh))
      
    } else if(identical(slope,"increasing")){
      g0 <- function(x){as.numeric(g(x, matrix(coef(object), 1, length(coef(object)))) - 
                                     g(0, matrix(coef(object), 1, length(coef(object)))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " + as.matrix(cbind(0, -1, 0, 0))[, notFixed]", sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " + as.matrix(cbind(0, -1, 0, 0, 0))[, notFixed]", sep="")
      }
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  if(def == "excess"){
    if(identical(slope,"decreasing")){
      g0 <- function(x){as.numeric((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(x, matrix(coef(object), 1, length(coef(object)))))/
                                     g(0, matrix(coef(object), 1, length(coef(object)))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(-1/(parmMat[,3]), 
                   -1/(1-parmMat[,3]), 
                   -1/(1-parmMat[,3]), -1/(1-parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 0, (", 
                     dpg, ") /(parmMat[,3]^2),
                   0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(-1/(parmMat[,3]), 
                   -1/(1-parmMat[,3]), -1/(1-parmMat[,3]), 
                   -1/(1-parmMat[,3]), -1/(1-parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 0, (", 
                     dpg, ") /(parmMat[,3]^2),
                   0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    } else if(identical(slope,"increasing")){
      g0 <- function(x){as.numeric((g(x, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(0, matrix(coef(object), 1, length(coef(object)))))/
                                     (1 - g(0, matrix(coef(object), 1, length(coef(object))))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(1-parmMat[,2]), 
                   1/(1-parmMat[,2]), 
                   1/(1-parmMat[,2]), 1/(1-parmMat[,2])))[, notFixed]",
                     " + as.matrix(cbind(0, (", 
                     dpg, " - 1)/(1-parmMat[,2])^2,
                   0, 0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(1-parmMat[,2]), 
                   1/(1-parmMat[,2]), 
                   1/(1-parmMat[,2]), 1/(1-parmMat[,2]), 1/(1-parmMat[,2])))[, notFixed]",
                     " + as.matrix(cbind(0, (", 
                     dpg, " - 1)/(1-parmMat[,2])^2,
                   0, 0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  if(def == "relative"){
    if(identical(slope,"decreasing")){
      g0 <- function(x){as.numeric((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(x, matrix(coef(object), 1, length(coef(object)))))/
                                     g(0, matrix(coef(object), 1, length(coef(object)))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]", dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(-1/(parmMat[,3]), 
                   -1/(parmMat[,3]), 
                   -1/(parmMat[,3]), -1/(parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 0, (", 
                     dpg, ") /(parmMat[,3]^2),
                   0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(-1/(parmMat[,3]), 
                   -1/(parmMat[,3]), -1/(parmMat[,3]), 
                   -1/(parmMat[,3]), -1/(parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 0, (", 
                     dpg, ") /(parmMat[,3]^2),
                   0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
      
    } else if(identical(slope,"increasing")){
      g0 <- function(x){as.numeric((g(x, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(0, matrix(coef(object), 1, length(coef(object)))))/
                                     (g(0, matrix(coef(object), 1, length(coef(object))))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,2]), 
                   1/(parmMat[,2]), 
                   1/(parmMat[,2]), 1/(parmMat[,2])))[, notFixed]",
                     " - as.matrix(cbind(0, (", 
                     dpg, " )/(parmMat[,2]^2),
                   0, 0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,2]), 
                   1/(parmMat[,2]), 
                   1/(parmMat[,2]), 1/(parmMat[,2]), 1/(parmMat[,2])))[, notFixed]",
                     " - as.matrix(cbind(0, (", 
                     dpg, ")/(parmMat[,2]^2),
                   0, 0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  if(def == "point"){
    stop(paste("Inverse regression not possible for def=point", sep=""))
  }
  
  if(def == "extra"){
    if(identical(slope,"decreasing")){
      g0 <- function(x){as.numeric((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(x, matrix(coef(object), 1, length(coef(object)))))/
                                     (g(0, matrix(coef(object), 1, length(coef(object)))) 
                                      - g(Inf, matrix(coef(object), 1, length(coef(object))))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,2]-parmMat[,3]), 
                   1/(parmMat[,2]-parmMat[,3]), 
                   1/(parmMat[,2]-parmMat[,3]), 1/(parmMat[,2]-parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (", dpg, " - parmMat[,2])/(parmMat[,2]-parmMat[,3])^2,
                     ( parmMat[,3] - (", dpg, "))/(parmMat[,2]-parmMat[,3])^2,
                     0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,2]-parmMat[,3]), 
                   1/(parmMat[,2]-parmMat[,3]), 1/(parmMat[,2]-parmMat[,3]), 
                   1/(parmMat[,2]-parmMat[,3]), 1/(parmMat[,2]-parmMat[,3])))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (", dpg, " - parmMat[,2])/(parmMat[,2]-parmMat[,3])^2,
                     ( parmMat[,3] - (", dpg, "))/(parmMat[,2]-parmMat[,3])^2,
                     0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    } else if(identical(slope,"increasing")){
      g0 <- function(x){as.numeric((g(x, matrix(coef(object), 1, length(coef(object)))) - 
                                      g(0, matrix(coef(object), 1, length(coef(object)))))/
                                     (g(Inf, matrix(coef(object), 1, length(coef(object)))) 
                                      - g(0, matrix(coef(object), 1, length(coef(object))))))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,3]-parmMat[,2]), 
                   1/(parmMat[,3]-parmMat[,2]), 
                   1/(parmMat[,3]-parmMat[,2]), 1/(parmMat[,3]-parmMat[,2])))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (", dpg, " - parmMat[,3])/(parmMat[,3]-parmMat[,2])^2,
                     ( parmMat[,2] - (", dpg, "))/(parmMat[,3]-parmMat[,2])^2,
                     0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind(1/(parmMat[,3]-parmMat[,2]), 
                   1/(parmMat[,3]-parmMat[,2]), 
                   1/(parmMat[,3]-parmMat[,2]), 1/(parmMat[,3]-parmMat[,2]), 1/(parmMat[,3]-parmMat[,2])))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (", dpg, " - parmMat[,3])/(parmMat[,3]-parmMat[,2])^2,
                     ( parmMat[,2] - (", dpg, "))/(parmMat[,3]-parmMat[,2])^2,
                     0, 0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  ###
  if((def == "hybridAdd" | def == "hybridExc") & !(backgType %in% c("hybridSD","absolute"))){
    stop(paste("When def = ", def, "backgType needs to be either hybridSD or absolute"))
  }
  
  if(def == "hybridAdd" & identical(backgType,"hybridSD")){
    if(identical(slope,"decreasing")){
      NumSD <- ifelse(is.na(backg), 2, backg)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric(pnorm((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                            NumSD*useSD-g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - 
                                     p0)}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," - (",dpg,"))/(",useSD,")))/",useSD,", 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,", 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    } else if(identical(slope,"increasing")){
      NumSD <- ifelse(is.na(backg), 2, backg)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric(1 - pnorm((NumSD*useSD + g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                                g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - p0)}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]", dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,", 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,", 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  if(def == "hybridAdd" & identical(backgType,"absolute")){
    if(is.na(backg)){
      stop(paste("No value of backg is provided", sep=""))
    }
    if(identical(slope,"decreasing")){
      NumSD <- abs((backg - g(0, matrix(coef(object), 1, length(coef(object)))))/useSD)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric(pnorm((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                            NumSD*useSD-g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - 
                                     p0)}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," - (",dpg,"))/(",useSD,")))/",useSD,", 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,", 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    } else if(identical(slope,"increasing")){
      NumSD <- abs((backg - g(0, matrix(coef(object), 1, length(coef(object)))))/useSD)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric(1 - pnorm((NumSD*useSD + g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                                g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - p0)}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]", dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,", 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,", 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,"))/",useSD,")))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/",useSD,",
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  
  ###
  
  if(def == "hybridExc" & identical(backgType,"hybridSD")){
    if(identical(slope,"decreasing")){
      NumSD <- ifelse(is.na(backg), 2, backg)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric((pnorm((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                             NumSD*useSD-g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - 
                                      p0)/(1-p0))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," - (",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/(",useSD,"*(1-",p0,")))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/(",useSD,"*(1-",p0,")))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
      
    } else if(identical(slope,"increasing")){
      NumSD <- ifelse(is.na(backg), 2, backg)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric((1 - pnorm((NumSD*useSD + g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                                 g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - p0)/(1-p0))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]", dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,"))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,"))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  if(def == "hybridExc" & identical(backgType,"absolute")){
    if(is.na(backg)){
      stop(paste("No value of backg is provided", sep=""))
    }
    if(identical(slope,"decreasing")){
      NumSD <- abs((backg - g(0, matrix(coef(object), 1, length(coef(object)))))/useSD)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric((pnorm((g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                             NumSD*useSD-g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - 
                                      p0)/(1-p0))}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]",dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," - (",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/(",useSD,"*(1-",p0,")))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (-dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,"))/(",useSD,"*(1-",p0,")))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (dnorm((parmMat[,3] - ",NumSD*useSD," -(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
      
    } else if(identical(slope,"increasing")){
      NumSD <- abs((backg - g(0, matrix(coef(object), 1, length(coef(object)))))/useSD)
      p0 <- 1-pnorm(NumSD)
      g0 <- function(x){as.numeric((1 - pnorm((NumSD*useSD + g(0, matrix(coef(object), 1, length(coef(object)))) - 
                                                 g(x, matrix(coef(object), 1, length(coef(object)))))/useSD) - p0))/(1-p0)}
      
      dpgh0 <- deparse(as.list(body(h))[[element.num]], width.cutoff = 500L)
      dpgh <- paste(dpgh0,collapse = "")
      dpg0 <- deparse(as.list(body(g))[[length(as.list(body(g)))]], width.cutoff = 500L)
      dpg <- gsub("cParm", "parmMat[,2]", dpg0)
      if(length(ParmVec0)==4){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,"))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0))[, notFixed]",
                     sep="")
      }
      if(length(ParmVec0)==5){
        pgh <- paste(dpgh, " * as.matrix(cbind((dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")), 
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                   (dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,"))))[, notFixed]",
                     " + as.matrix(cbind(0, 
                     (-dnorm((",NumSD*useSD,"+parmMat[,2]-(",dpg,"))/(",useSD,")))/(",useSD,"*(1-",p0,")),
                     0,0,0))[, notFixed]",
                     sep="")
      }
      
      body(h)[[element.num]] <- as.call(str2lang(pgh))
    }
  }
  
  ###  
  respType <- object[["type"]]
  if (identical(respType, "continuous"))
  {
    tquan <- qt(1 - (1 - level)/2, df.residual(object))   
  } else {
    tquan <- qnorm(1 - (1 - level)/2)
  }
  
  
  
  #  g1<-function(x){as.numeric(g(x, matrix(coef(object), 1, length(coef(object)))))}
  if(sandwich.vcov){
    j1<-function(x){sqrt(as.vector(h(x, matrix(coef(object), 1, length(coef(object)))))%*%
                           sandwich::sandwich(object)%*%
                           as.vector(h(x, matrix(coef(object), 1, length(coef(object))))))}
  } else{
    j1<-function(x){sqrt(as.vector(h(x, matrix(coef(object), 1, length(coef(object)))))%*%
                           vcov(object)%*%
                           as.vector(h(x, matrix(coef(object), 1, length(coef(object))))))}
  }
  my.fun1<-function(x){g0(x) + tquan * j1(x)}
  my.fun2<-function(x){g0(x) - tquan * j1(x)}
  
  objDL <- object[["dataList"]][["names"]]
  #colnames(newData0) <- c(objDL[["dName"]], objDL[["cName"]])
  maxdose <- extFactor * max(object[["dataList"]][["dose"]])
  
  
  rootFctBMD <- function(x) 
  {
    newData <- data.frame(x, catLev)
    colnames(newData) <- c(objDL[["dName"]], objDL[["cName"]])
    #    print(c(x, predict(object, newData, interval = intType, level = level)[2] - yval))
    g0(x) - bmr
  }
  
  rootFct1 <- function(x) 
  {
    newData <- data.frame(x, catLev)
    colnames(newData) <- c(objDL[["dName"]], objDL[["cName"]])
    my.fun1(x) - bmr
  }
  
  rootFct2 <- function(x) 
  {
    newData <- data.frame(x, catLev)
    colnames(newData) <- c(objDL[["dName"]], objDL[["cName"]])
    my.fun2(x) - bmr
  }
  
  ResMat <- matrix(NA,1,3)
  colnames(ResMat)<-c("BMD","BMDL","BMDU")
  BMDLevel<-try(uniroot(rootFctBMD, c(0, maxdose)), silent = TRUE)
  if(inherits(BMDLevel,"try-error")){
    stop(paste("No solution possible for this BMR", sep=""))
  } else {
    ResMat[1,1] <- BMDLevel$root
    
    upLevel<-try(uniroot(rootFct1, c(0, maxdose)), silent = TRUE)
    ResMat[1,2]<- ifelse(inherits(upLevel,"try-error"), 0, upLevel$root)
    
    lowLevel<-try(uniroot(rootFct2, c(0, maxdose)), silent = TRUE)
    ResMat[1,3]<- ifelse(inherits(lowLevel,"try-error"), Inf, lowLevel$root)
  }
  return(ResMat)
}

