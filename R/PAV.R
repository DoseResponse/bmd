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
