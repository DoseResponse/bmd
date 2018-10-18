aggregateBinomial <- function(object,data){
  df <- aggregate(cbind(data[,strsplit(as.character(object[[2]]),"/")[[2]]],
                        data[,strsplit(as.character(object[[2]]),"/")[[3]]]) ~
                    data[,as.character(object[[3]])],FUN = sum)
  colnames(df) <- c(as.character(object[[3]]),
                    strsplit(as.character(object[[2]]),"/")[[2]],
                    strsplit(as.character(object[[2]]),"/")[[3]])
  df
  
}
