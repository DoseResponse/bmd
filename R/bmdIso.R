#' Benchmark dose estimation from isotonic regression
#' 
#' The function monotonizes the sequence of response values based on the
#' pool-adjacent violators algorithm and based on these use linear
#' interpolating splines to build an isotonic regression model. From this model
#' a benchmark dose is estimated.
#' 
#' Formula should be specified as in the form number/total ~ dose for binomial
#' data. For details about the use of isotonic regression for BMD estimation
#' see Piegorsch et al. (20014) and Lin et al. (2015).
#' 
#' @param formula an object of class "formula" expressing the dose-response
#' relationship. Details of model specification are given under 'Details'
#' @param data data frame containing the variables in the model
#' @param type character string specifying the type of data used in the model.
#' "continuous", "binomial" or "Poisson"
#' @param bmr numeric value of benchmark response level for which to calculate
#' the benchmark dose
#' @param p0 background probability for hybrid definitions
#' @param backgType character string specifying how the background level is
#' specified. For binomial data the options are "modelBased" and "absolute".
#' For continuous data the options are "absolute", "hybridSD" and
#' "hybridPercentile"
#' 
#' "modelBased" - the background level is obtained from the model as the level
#' for dose 0: p0 = f(0)
#' 
#' "absolute" - the background level is specified by the user through the backg
#' argument: p0 = backg for binomial response and for the "relative" and
#' "added" definition for continuous response.  p0 = 1 - phi((back -
#' f(0))/sigma) for "hybridExc" and "hybridAdd" definitions.
#' 
#' "hybridSD" - the background risk is specified by the user in terms of number
#' of SDs from the mean of the control group.  p0 = 1 - phi(((backg*sigma +
#' f(0)) - f(0))/sigma) = 1 - phi(backg), where phi is the normal distribution
#' function and sigma is the SD for the control group.  "hybridPercentile" -
#' the background risk is specified by the user in terms of percentile from the
#' control group distribution (assuming a normal distribution).  p0 = 1 -
#' phi((x0 - f(0))/sigma) = 1 - backg.  where x0 is the level for which the
#' response is considered adverse, phi is the normal distribution function and
#' sigma is the SD for the control group
#' @param backg numeric value specifying the background level. Defaults to 0
#' for "absolute" background risk for binomial response (1 for decreasing
#' dose-response models), 2 SD for "hybridSD" background and 0.9 for
#' "hybridpercentile"
#' @param def character string specifying the definition of the benchmark dose
#' to use in the calculations. "excess", "additional" and "point" are for
#' binomial response whereas "relative", "added", "hybridExc" (excess hybrid),
#' "hybridAdd" (additional hybrid), and "point" are for continuous response.
#' "relative", "extra", and "point" are for count response data.
#' 
#' "excess" - BMR is defined as: BMR = (f(BMD) - p0)/(1 - p0).  Works for
#' binomial response. BMR should be between 0 and 1.
#' 
#' "additional" - BMR is defined as: BMR = f(BMD) - p0.  Works for binomial
#' response. BMR should be between 0 and 1.
#' 
#' "point" - The response level for which to find BMD is directly defined
#' through the BMR level: BMR = f(BMD). Works for binomial, count and
#' continuous response
#' 
#' "relative" - BMR is defined as: BMR = (f(BMD) - p0)/p0.  Works for count and
#' continuous response
#' 
#' "added" - BMR is defined as: BMR= f(BMD) + p0.  Works for continuous
#' response
#' 
#' "hybridExc" - BMR is defined as: BMR = (1 - phi((x0 - f(BMD))/sigma) - p0)/
#' (1- p0), where x0 is the level for which the response is considered adverse,
#' phi is the normal distribution function and sigma is the SD for the control
#' group.  Works for continuous response
#' 
#' "hybridAdd" - BMR is defined as: BMR = 1 - phi((x0 - f(BMD))/sigma) - p0,
#' where x0 is the level for which the response is considered adverse, phi is
#' the normal distribution function and sigma is the SD for the control group.
#' Works for continuous response
#' @param display logical. If TRUE the results are displayed; otherwise they
#' are not
#' @return A matrix with two columns, one containing BMD and the other
#' containing BMDL.
#' @author Signe M. Jensen
#' @references Piegorsch, W. W., Xiong, H., Bhattacharya, R. N., & Lin, L.
#' (2014). Benchmark dose analysis via nonparametric regression modeling. Risk
#' Analysis, 34(1), 135-151
#' 
#' Lin, L., Piegorsch, W. W. and Bhattacharya R. (2015). Nonparametric
#' benchmark dose estimation with continuous dose-response data. Scandinavian
#' Journal of Statistics, 42, 713-731
#' @keywords nonparametric isotonic regression
#' @examples
#' 
#' ## Data on tumor incidence in rats after exposure to formaldehyde, from Piegorsch et al. (2014)
#' formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
#'                           tumor.incidence = c(0, 0, 0, 3, 21, 150),
#'                           total = c(122, 27, 126, 113, 34, 182))
#'                           
#' # Estimating BMD from isotonic regression using excess risk definition and a BMR=0.1
#' bmdIso(tumor.incidence/total ~ conc, 
#'       data=formaldehyde, 
#'       type="binomial",
#'       bmr=0.1,
#'       backgType = "modelBased",
#'       def = "excess")
#'       
#'       
#' ## Data on root length in ryegrass after exposure to ferulic acid
#' require(drc)
#' require(drcData)
#' data(ryegrass)
#' 
#' # As isotonic regression only wors for increasing dose-response relationship 
#' # the association is turned
#' ryegrass1<-ryegrass
#' ryegrass1$rootl<-100-ryegrass1$rootl
#' 
#' # Estimating BMD from isotonic regression using relative risk definition
#' # and a BMR=0.05
#' bmdIso(rootl ~ conc, 
#'       data=ryegrass1, 
#'       type="continuous",
#'       bmr=0.05,
#'       backgType = "modelBased",
#'       def = "relative")
#' @export
bmdIso <- function(formula, data, type, bmr, p0, backgType = c("modelBased", "absolute","hybridSD","hybridPercentile"), backg=NA, def = c("excess", "additional", "relative", "added", "hybridExc", "hybridAdd", "point"), display=FALSE){
  object <- formula
  PAV.p <- PAV(object, data, type)
  n <- as.numeric(table(data[, paste(object[[3]])]))
  if(type=="continuous"){
    #sigma.sq <- sum((data[,paste(object[[2]])]-rep(PAV.p,n))^2)/length(data[,paste(object[[2]])])
    sigma.sq <- sd(data[data[,paste(object[[3]])]==0,paste(object[[2]])])
  }
  Dose <- sort(unique(data[,paste(object[[3]])]))
  f0 <- min(PAV.p) 
  if (missing(def)) {
    stop(paste("def is missing", sep=""))
  }
  if (!(def %in% c("excess", "additional", "relative", "added", "hybridExc", "hybridAdd", "point"))) {
    stop(paste("Could not recognize def", sep=""))
  }
  if (identical(backgType,"modelBased")) {
    background <- f0
  } else if (identical(backgType,"absolute") & !(def %in% c("relative"))) {
    background <- backg
  } else if (identical(backgType,"absolute") & !(def %in% c("hybridExc","hybridAdd"))) {
    background <- ifelse(is.na(backg),0,backg)
  } else if (identical(backgType,"hybridSD")) {
    background <- ifelse(is.na(backg), 1-pnorm(2), 1-pnorm(backg))
  } else if (identical(backgType,"absolute") & 
             (identical(def,"hybridExc") | identical(def,"hybridAdd") )) {
    background <- ifelse(is.na(backg), 
                         1-pnorm(2),
                         1-pnorm((backg-f0)/sigma.sq))
  } else {
    background <- ifelse(is.na(backg),1-0.9,1-backg)
  }
  def <- match.arg(def)
  if (identical(type, "binomial")) {
    Cq <- switch(def, 
                 excess = bmr * (1 - background) + background, 
                 additional = bmr + background, 
                 point = bmr)
  }
  if (identical(type, "binomial") & (def %in% c("relative","added", "hybridExc","hybridAdd"))) {
    stop(paste("\"",def, "\" is not available for quantal data", sep=""))
  }
  if (type %in% c("Poisson","negbin1","negbin2")) {
      Cq <- switch(def,
                   relative = bmr * background + background,
                   extra = bmr*abs(diff(predict(object, data.frame(c(0, Inf))))) + background,
                   point = bmr)
    }
  if (type %in% c("Poisson","negbin1","negbin2") & (def %in% c("excess","additional","added","hybridExc","hybridAdd"))) {
    stop(paste("\"",def, "\" is not available for count data", sep=""))
  }
  if (identical(type, "continuous") & (def %in% c("excess", "additional"))) {
    stop(paste("\"",def, "\" is not available for continuous data", sep=""))
  }
  if (identical(type, "continuous")) {
    Cq <- switch(def, relative = bmr * background + background, 
                 added = bmr + background,
                 point = bmr,
                 hybridAdd = sigma.sq * 
                   (qnorm(1 - background) - qnorm(1 - (background + bmr))) + 
                   f0,
                 hybridExc = sigma.sq * 
                   (qnorm(1 - background) - qnorm(1 + background - (1 - background)*bmr)) + 
                   f0)
  } 
  if(sum(Cq == PAV.p)>0){
    BMD <- max(Dose[Cq == PAV.p])
  } else if( length(PAV.p) >  sum(as.numeric(Cq > PAV.p)) ){
    crit <- sum(as.numeric(Cq > PAV.p))
    BMD <- Dose[crit] + (Cq-PAV.p[crit])*(Dose[crit+1]-Dose[crit])/(PAV.p[crit+1]-PAV.p[crit])  
  } else if(sum(as.numeric(Cq < PAV.p))==0){
    BMD <- 0
  } else if(sum(as.numeric(Cq < PAV.p)) == length(PAV.p)){
    BMD <- Dose[length(PAV.p)]
  }
  if (display) {
    cat("Effective response level: ", Cq )
  }
  BMD
}
