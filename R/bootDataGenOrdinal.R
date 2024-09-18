bootDataGenOrdinal <- function(object, R = 500, bootType = c("nonparametric", "parametric", "model", "hierarchical")){
  bootType <- match.arg(bootType)
  
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
      for(j in 1:length(object$levels)){
        if(!(object$levels[j] %in% colnames(df))){
          df[,object$levels[j]]<-0
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
        prop0[prop0==0] <- (1/length(object$levels)^2)/(data.size+1/length(object$levels)) # (1/4)/(data.size + 1/2)
        prop0[prop0==1] <- (data.size + 1/length(object$levels)^2)/(data.size+1/length(object$levels)) # (data.size + 1/4)/(data.size+1/2)
        sampled.expand$variable[sampled.expand$row.orig==j] <-
          #rep(unlist(object$levels), as.numeric(rmultinom(1, size = data.size, prob = prop0)))
          unlist(sample(object$levels, size = data.size, replace = TRUE, prob = prop0))
      }
      columns.rem <- colnames(subset(sampled.expand, select=-c(variable,value,row.orig)))
      df <- reshape2:::dcast(sampled.expand, as.formula(paste(paste(columns.rem, collapse = "+")," ~ variable")), length)
      for(j in 1:length(object$levels)){
        if(!(object$levels[j] %in% colnames(df))){
          df[,object$levels[j]]<-0
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
      sampled.expand[, "variable"] <- sapply(data.e[,object$dose], function(x) unlist(sample(object$levels, size = 1, replace = TRUE, object$pFun(x))))
      columns.rem <- colnames(subset(sampled.expand, select=-c(variable,value,row.orig)))
      df <- reshape2:::dcast(sampled.expand, as.formula(paste(paste(columns.rem, collapse = "+")," ~ variable")), length)
      for(j in 1:length(object$levels)){
        if(!(object$levels[j] %in% colnames(df))){
          df[,object$levels[j]]<-0
        }
      }
      tmp.data[[i]] <- df
    }
  }
  if(bootType == "hierarchical"){
    if(is.null(object$blocks)){
      stop('ordinal dose-response model does not include blocks. Hierarchical resampling is not possible.')
      # cat("\n", 'Argument "block" needs to be specified for hierarchical resampling.', "\n")
      # tmp.data <- NULL
    } 
    
    resample_fun <- function(levels, dose, weights, blocks, data){
      data %>%
        dplyr::mutate(row.orig = 1:n()) %>% 
        dplyr::group_by(.data[[dose]]) %>% 
        dplyr::slice_sample(prop = 1, weight_by = eval(parse(text=weights)), replace = TRUE) %>% 
        tidyr::pivot_longer(levels) %>% 
        dplyr::group_by(row.orig) %>% 
        tidyr::uncount(value) %>% 
        dplyr::group_by(row.orig) %>% 
        dplyr::slice_sample(prop = 1, replace = TRUE) %>% 
        dplyr::group_by(.data[[dose]], .data[[weights]], .data[[blocks]], row.orig) %>% 
        dplyr::count(name) %>% 
        dplyr::mutate(total = sum(n)) %>% 
        tidyr::pivot_wider(names_from = name, values_from = n, values_fill = 0)
    }
    
    tmp.data <- list()
    for (i in 1:R){
      tmp.data[[i]] <- resample_fun(object$levels, object$dose, object$weights, object$blocks, object$data)
    }
  }
  tmp.data
}
