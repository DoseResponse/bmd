expandBinomial <- function(data, number, total, dose){
  df<-data.frame(number = c(rep(rep(1,length(data[,number])), data[,number]),
                            rep(rep(0,length(data[,number])), 
                                data[,total]-data[,number])),
                 dose = c(rep(data[,dose], data[,number]),
                          rep(data[,dose], data[,total]-data[,number])),
                 total = 1)
  colnames(df)<-c(number, dose, total)
  df
}
