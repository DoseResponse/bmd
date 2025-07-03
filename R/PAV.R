#' Pool-adjacent-violators monotonizing
#' 
#' The function monotonizes a sequence of probabilities or means based on the
#' pool-adjacent-violators algorithm.
#' 
#' For details on how the pool-adjacent-violators algorithm is defined see
#' Silvapulle and Sen (2004).
#' 
#' Formula should be specified as in the form number/total ~ dose for binomial
#' data and response ~ for continuous data.
#' 
#' @param formula an object of class "formula" expressing dose-response
#' relationship. Details of model specification are given under 'Details'
#' @param data data frame containing the variables in the formula
#' @param type character string specifying the type of data used in the model,
#' "continuous" or "binomial" or "Poisson"
#' @return A vector containing the monotonized sequence.
#' @author Signe M. Jensen
#' @references Silvapulle, M. J. and Sen, P. K. (2004). Constrained statistical
#' inference: order, inequality, and shape constraints. New York: John Wiley &
#' Sons.
#' @keywords nonparametric isotonic regression
PAV<-function(formula,data,type){
  object <- formula
  if( identical(type,"binomial")){
    N <- length( data[, paste(object[[3]]) ])
    Events <- data[, strsplit(as.character(object[[2]]),"/")[[2]] ]
    Total <- data[, strsplit(as.character(object[[2]]),"/")[[3]] ]
    Dose <- data[,paste(object[[3]])]
    PAV.p <- rep(NA,N)
    for(i in 1:N){
      tmp2 <- rep(NA,length(1:i))
      for(u in 1:i){
        tmp1 <- rep(NA,length(1:N))
        for(v in i:N){
          tmp1[v] <- sum(Events[u:v])/sum(Total[u:v])
        }
        tmp2[u] <- min(tmp1,na.rm = TRUE)
      }
      PAV.p[i] <- max(tmp2)
    }
  }
  if(type %in% c("continuous","Poisson","negbin1","negbin2")){
    N <- length( unique(data[, paste(object[[3]]) ]))
    Response.m<-aggregate(data[,as.character(object[[2]])] ~data[,as.character(object[[3]])],FUN=mean)[,2]
    n <- as.numeric(table(data[, paste(object[[3]])]))
    PAV.p<-rep(NA,N)
    for(i in 1:N){
      tmp2<-rep(NA,length(1:i))
      for(u in 1:i){
        tmp1<-rep(NA,length(1:N))
        for(v in i:N){
          tmp1[v]<-sum(n[u:v]*Response.m[u:v])/sum(n[u:v])
        }
        tmp2[u]<-min(tmp1,na.rm = TRUE)
      }
      PAV.p[i]<-max(tmp2)
    }
  }
  PAV.p
}
