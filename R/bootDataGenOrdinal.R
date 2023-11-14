bootDataGenOrdinal <- function(object, R = 500, bootType = "nonparametric"){
  if(!(bootType %in% c("nonparametric", "parametric", "model"))) cat("\n", "Error: bootType not defined", "\n")
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
      df <- reshape2:::dcast(sampled.expand, as.formula(paste(paste(columns.rem, collapse = "+")," ~ variable")), length)
      for(j in 1:length(object$categories)){
        if(!(object$categories[j] %in% colnames(df))){
          df[,object$categories[j]]<-0
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
        prop0[prop0==0] <- (1/length(object$categories)^2)/(data.size+1/length(object$categories)) # (1/4)/(data.size + 1/2)
        prop0[prop0==1] <- (data.size + 1/length(object$categories)^2)/(data.size+1/length(object$categories)) # (data.size + 1/4)/(data.size+1/2)
        sampled.expand$variable[sampled.expand$row.orig==j] <-
          #rep(unlist(object$categories), as.numeric(rmultinom(1, size = data.size, prob = prop0)))
          unlist(sample(object$categories, size = data.size, replace = TRUE, prob = prop0))
      }
      columns.rem <- colnames(subset(sampled.expand, select=-c(variable,value,row.orig)))
      df <- reshape2:::dcast(sampled.expand, as.formula(paste(paste(columns.rem, collapse = "+")," ~ variable")), length)
      for(j in 1:length(object$categories)){
        if(!(object$categories[j] %in% colnames(df))){
          df[,object$categories[j]]<-0
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
      sampled.expand[, "variable"] <- sapply(data.e[,object$dose], function(x) unlist(sample(object$categories, size = 1, replace = TRUE, object$pFun(x))))
      columns.rem <- colnames(subset(sampled.expand, select=-c(variable,value,row.orig)))
      df <- reshape2:::dcast(sampled.expand, as.formula(paste(paste(columns.rem, collapse = "+")," ~ variable")), length)
      for(j in 1:length(object$categories)){
        if(!(object$categories[j] %in% colnames(df))){
          df[,object$categories[j]]<-0
        }
      }
      tmp.data[[i]] <- df
    }
  }
  tmp.data
}
