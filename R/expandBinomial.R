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
