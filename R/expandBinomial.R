#' Utility function
#' 
#' Utility function for expanding binomial data
#' 
#' The aim to provide an R package calculating the benchmark dose (BMD) and the
#' lower limit of the corresponding 95\% confidence interval (BMDL) for
#' continuous and quantal dose-response data for a range of dose-response
#' models based on the available definitions of the benchmark dose concepts.
#' 
#' @param data a data.frame
#' @param number the name of the column in the data set containing the number
#' of affected individuals per dose level
#' @param total the name of the column in the data set containing the total
#' number of individuals per dose level
#' @param dose the name of the column in the data set containing the dose
#' values
#' @param curveid the name of the column in the data set specifying the curveid
#' (if available)
#' @return data.frame
#' @author Signe M. Jensen and Jens Riis Baalkilde
#' @keywords models nonlinear
#' @examples \dontrun{
#' }
#'  
#' 
expandBinomial <- function(data, number, total, dose, curveid = character(0)){
  if(length(curveid) == 0){
    df<-data.frame(number = c(rep(rep(1,length(data[,number])), data[,number]),
                              rep(rep(0,length(data[,number])), 
                                  data[,total]-data[,number])),
                   dose = c(rep(data[,dose], data[,number]),
                            rep(data[,dose], data[,total]-data[,number])),
                   total = 1)
    colnames(df)<-c(number, dose, total)
  } else {
    df<-data.frame(number = c(rep(rep(1,length(data[,number])), data[,number]),
                              rep(rep(0,length(data[,number])), 
                                  data[,total]-data[,number])),
                   dose = c(rep(data[,dose], data[,number]),
                            rep(data[,dose], data[,total]-data[,number])),
                   total = 1,
                   curveid = c(rep(data[,paste0("orig.", curveid)], data[,number]),
                               rep(data[,paste0("orig.", curveid)], data[,total]-data[,number])))
    colnames(df)<-c(number, dose, total, curveid)
  }
  df
}
