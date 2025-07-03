#' Benchmark dose estimation from isotonic regression
#' 
#' The function estimates BMD and BMDL using bootstrap based on isotonic
#' regression.
#' 
#' BMD and BMDL is defined as the median and the 5th percentile in the
#' bootstrap distribution, respectively.
#' 
#' Formula should be specified as in the form number/total ~ dose for binomial
#' data.
#' 
#' Bootstrapping with the argument boot = "resample" is done by sampling with
#' replacement from the original data set. Bootstrapping with the argument boot
#' = "pseudorandom" is done by sampling from norm(mean(Y_i),sd(Y_0)), assuming
#' equal variance between groups, in case of continuous data. For binomial
#' data, each bootstrap data set is sampled from binom(N_i,Y_i/N_i). In case of
#' Y_i = 0 or Y_i = N_i shrinkage is used to avoid that the resampling always
#' produces 0 or 1, respectively. In this case data is sampled from
#' binom(N_i,(Y_i+1/3)/(N_i+1/3)).
#' 
#' All sampling is made within dose groups.
#' 
#' For details about the use of isotonic regression for BMD estimation see
#' Piegorsch et al. (20014) and Lin et al. (2015).
#' 
#' @param formula an object of class "formula" expressing dose-response
#' relationship. Details of model specification are given under 'Details'
#' @param data data frame containing the variables in the model
#' @param type character string specifying the type of data used in the model.
#' "continuous" or "binomial" or "Poisson"
#' @param bmr numeric value of benchmark response level for which to calculate
#' the benchmark dose
#' @param R number of bootstrap samples
#' @param boot character string specifying type of bootstrap used. "resample"
#' or "pseudorandom". Only option for count data is "resample"
#' @param backgType character string specifying how the background level is
#' specified. For binomial data the options are "modelBased" and "absolute".
#' For continuous data the options are "modelBased", "absolute", "hybridSD" and
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
#' @keywords nonparametric isotonic regression bootstrap
#' @examples
#' 
#' ## Data on tumor incidence in rats after exposure to formaldehyde, from Piegorsch et al. (2014)
#' formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
#'                           tumor.incidence = c(0, 0, 0, 3, 21, 150),
#'                           total = c(122, 27, 126, 113, 34, 182))
#'                           
#' # BMD and BMDL from isotonic regression using excess risk definition and a BMR=0.1
#' bmdIsoBoot(tumor.incidence/total ~ conc, 
#'       data=formaldehyde, 
#'       type="binomial",
#'       bmr=0.1,
#'       backgType = "modelBased",
#'       def = "excess")
#'       
#'       
#' ## Data on root length in ryegrass after exposure to ferulic acid
#' require(drc)
#' data(ryegrass)
#' 
#' # As isotonic regression only wors for increasing dose-response relationship
#' # the association is turned
#' ryegrass1<-ryegrass
#' ryegrass1$rootl<-100-ryegrass1$rootl
#' 
#' # Estimating BMD from isotonic regression using relative risk definition and a BMR=0.05
#' bmdIsoBoot(rootl ~ conc, 
#'       data=ryegrass1, 
#'       type="continuous",
#'       bmr=0.05,
#'       backgType = "modelBased",
#'       def = "relative", R = 100)
#' 
bmdIsoBoot <- function(formula, data, type, bmr, R=1000, boot="resample", 
                       backgType = c("modelBased", "absolute","hybridSD","hybridPercentile"), 
                       backg=NA, 
                       def = c("excess", "additional", "relative", "added", "hybridExc", "hybridAdd", "point")){
  object <- formula
  if (type %in% c("Poisson","negbin1","negbin2") & boot!="resample") {
    stop(paste("\"",type,"\" only works with resample bootstrap \"", sep=""))
  }
  if(boot=="resample"){
    if(type=="binomial"){
      data.e<-expandBinomial(data, 
                          number = strsplit(as.character(object[[2]]),"/")[[2]],
                          total = strsplit(as.character(object[[2]]),"/")[[3]],
                          dose = as.character(object[[3]]))
  data.e[,"row.num"]<-1:dim(data.e)[1]
  tmp.data <- list()
  for(i in 1:R){
    sampled.expand <- data.e[as.numeric(unlist(aggregate(row.num ~ data.e[,as.character(object[[3]])], data=data.e, 
                                                        FUN=function(x) sample(x,replace=TRUE))[[2]])),]
    tmp.data[[i]] <- aggregateBinomial(object, sampled.expand)
  }
    }
    if(type %in% c("continuous","Poisson","negbin1","negbin2")){
      data.e<-data
      data.e[,"row.num"]<-1:dim(data.e)[1]
      data.e[,"dose"]<-data.e[,as.character(object[[3]])]
      tmp.data <- list()
      for(i in 1:R){
        tmp.data[[i]] <- data.e[as.numeric(unlist(aggregate(row.num ~ dose, data=data.e, 
                                                             FUN=function(x) sample(x,replace=TRUE))[[2]])),]
         }
      }
  } else if(boot=="pseudorandom"){
    if(type=="binomial"){
    Y <- data[,strsplit(as.character(object[[2]]),"/")[[2]]]
    N <- data[,strsplit(as.character(object[[2]]),"/")[[3]]]
    shrinks <- which(Y==N | Y==0)
    Y[shrinks] <- Y[shrinks]+0.25
    N[shrinks] <- N[shrinks]+0.5
    prob <- rep(Y/N,N)
    tmp.data <- list()
    for(i in 1:R){
      sampled.expand <- data.frame(number = rbinom(length(prob),1,prob), 
                                   dose = rep(data[,as.character(object[[3]])],N), 
                                   total = 1)
      df <- aggregateBinomial(number/total~dose, sampled.expand)
      colnames(df) <- c(as.character(object[[3]]),
                                        strsplit(as.character(object[[2]]),"/")[[2]],
                                        strsplit(as.character(object[[2]]),"/")[[3]])
      tmp.data[[i]] <- df  
    }
    }
    if(type=="continuous"){
      mean.Y <- aggregate(data[,as.character(object[[2]])]~data[,as.character(object[[3]])],FUN=function(x) mean(x,na.rm=TRUE))[,2]
      sd.Y <- aggregate(data[,as.character(object[[2]])]~data[,as.character(object[[3]])],FUN=function(x) sd(x,na.rm=TRUE))[,2]
      Dose<- aggregate(data[,as.character(object[[2]])]~data[,as.character(object[[3]])],FUN=function(x) length(!is.na(x)))[,1]
      N.dose<- aggregate(data[,as.character(object[[2]])]~data[,as.character(object[[3]])],FUN=function(x) length(!is.na(x)))[,2]
      tmp.data <- list()
      for(i in 1:R){
        sampled <- data.frame(y = rnorm(sum(N.dose),mean=rep(mean.Y,N.dose),sd=rep(sd.Y,N.dose)), 
                                     dose = rep(Dose,N.dose))
        colnames(sampled) <- c(as.character(object[[2]]), as.character(object[[3]]))
        tmp.data[[i]] <- sampled
        }
    }
  }
    bmd.list <- lapply(tmp.data, function(x){
    bmdIso(object, data=x, type=type, bmr = bmr, backgType = backgType, backg=backg,def=def)})
  
    resMat <- matrix(NA,1,2)
    resMat[1,1] <- bmdIso(object, data=data, type=type, bmr = bmr, backgType = backgType, backg=backg,def=def) # quantile(unlist(bmd.list),0.5)
    resMat[1,2] <- quantile(unlist(bmd.list),0.05)
    colnames(resMat) <- c("BMD", "BMDL")
    rownames(resMat) <- c("")
    cat("\n\n")
    resMat    
}

