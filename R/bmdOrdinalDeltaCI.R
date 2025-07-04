# bmdOrdinalDeltaCI <- function(object, bmr, backgType, backg, def, level){
#   if(substr(object$drmList[[1]]$fct$name, 1,2) %in% c("FP") || (object$drmList[[1]]$fct$name %in% c("LN.3", "LN.4"))){
#     warning("delta CI not available for FPL or LN.3 and LN.4 models.")
#     CI <- NA
#   } else {
#     if(backgType == "modelBased"){
#       # log-logistic
#       if(object$drmList[[1]]$fct$name == "LL.2"){
#         closedExpression <- paste0("exp(log(1/",bmr,"-1)/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "LL.3"){
#         closedExpression <- paste0("exp(log(d/",bmr,"-1)/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "LL.3u"){
#         if(def == "excess"){
#           closedExpression <- paste0("exp(log(1/",bmr,"-1)/b)*e")
#         } else if(def == "additional"){
#           closedExpression <- paste0("exp(log((1-c)/",bmr,"-1)/b)*e")
#         } else if(def == "point"){
#           closedExpression <- paste0("exp(log((1-c)/(",bmr,"-c)-1)/b)*e")
#         }
#       } else if(object$drmList[[1]]$fct$name == "LL.4"){
#         if(def == "excess"){
#           closedExpression <- paste0("exp(log((d-c)/((1-c)*",bmr,")-1)/b)*e")
#         } else if(def == "additional"){
#           closedExpression <- paste0("exp(log((d-c)/",bmr,"-1)/b)*e")
#         } else if(def == "point"){
#           closedExpression <- paste0("exp(log((d-c)/(",bmr,"-c)-1)/b)*e")
#         }
#       }
#       # log-normal
#       else if(object$drmList[[1]]$fct$name == "LN.2"){
#         constant <- qnorm(bmr)
#         closedExpression <- paste0("exp(",constant,"/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "LN.3"){
#         constant <- qnorm(bmr)
#         closedExpression <- paste0("exp(",constant,"/b)*e") # this is wrong
#       } else if(object$drmList[[1]]$fct$name == "LN.3u"){
#         constant <- qnorm(bmr)
#         closedExpression <- paste0("exp(",constant,"/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "LN.4"){
#         constant <- qnorm(bmr)
#         closedExpression <- paste0("exp(",constant,"/b)*e") # this is wrong
#       }
#       # weibull 1
#       else if(object$drmList[[1]]$fct$name == "W1.2"){
#         closedExpression <- paste0("(-log(",bmr,"))^(1/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "W1.3"){
#         closedExpression <- paste0("(-log(",bmr,"/d))^(1/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "W1.3u"){
#         if(def == "excess"){
#           closedExpression <- paste0("(-log(",bmr,"))^(1/b)*e")
#         } else if(def == "additional"){
#           closedExpression <- paste0("(-log(",bmr,"/(1-c)))^(1/b)*e")
#         } else if(def == "point"){
#           closedExpression <- paste0("(-log((",bmr,"-c)/(1-c)))^(1/b)*e")
#         }
#       } else if(object$drmList[[1]]$fct$name == "W1.4"){
#         if(def == "excess"){
#           closedExpression <- paste0("(-log((1-c)*",bmr,"/(d-c)))^(1/b)*e")
#         } else if(def == "additional"){
#           closedExpression <- paste0("(-log(",bmr,"/(d-c)))^(1/b)*e")
#         } else if(def == "point"){
#           closedExpression <- paste0("(-log((",bmr,"-c)/(d-c)))^(1/b)*e")
#         }
#       } 
#       # weibull 2
#       else if(object$drmList[[1]]$fct$name == "W2.2"){
#         closedExpression <- paste0("(-log(1-",bmr,"))^(1/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "W2.3"){
#         closedExpression <- paste0("(-log(1-",bmr,"/d))^(1/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "W2.3u"){
#         if(def == "excess"){
#           closedExpression <- paste0("(-log(1-",bmr,"))^(1/b)*e")
#         } else if(def == "additional"){
#           closedExpression <- paste0("(-log(1-",bmr,"/(1-c)))^(1/b)*e")
#         } else if(def == "point"){
#           closedExpression <- paste0("(-log(1-(",bmr,"-c)/(1-c)))^(1/b)*e")
#         }
#       } else if(object$drmList[[1]]$fct$name == "W2.4"){
#         if(def == "excess"){
#           closedExpression <- paste0("(-log(1-(1-c)*",bmr,"/(d-c)))^(1/b)*e")
#         } else if(def == "additional"){
#           closedExpression <- paste0("(-log(1-",bmr,"/(d-c)))^(1/b)*e")
#         } else if(def == "point"){
#           closedExpression <- paste0("(-log(1-(",bmr,"-c)/(d-c)))^(1/b)*e")
#         }
#       } 
#     } else if(backgType == "absolute"){
#       p0 <- ifelse(is.na(backg), 0, backg)
#       if(def == "excess"){
#         z <- (1-p0)*bmr + p0
#       } else if(def == "additional"){
#         z <- bmr + p0
#       }
#       # log-logistic
#       if(object$drmList[[1]]$fct$name == "LL.2"){
#         closedExpression <- paste0("exp(log(1/",z, "-1)/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "LL.3"){
#         closedExpression <- paste0("exp(log(d/",z,"-1)/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "LL.3u"){
#         closedExpression <- paste0("exp(log((1-c)/(",z,"-c)-1)/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "LL.4"){
#         closedExpression <- paste0("exp(log((d-c)/(",z,"-c)-1)/b)*e")
#       }
#       # log-normal
#       else if(object$drmList[[1]]$fct$name == "LN.2"){
#         constant <- qnorm(z)
#         closedExpression <- paste0("exp(",constant,"/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "LN.3"){
#         constant <- qnorm(bmr)
#         closedExpression <- paste0("exp(",constant,"/b)*e") # this is wrong
#       } else if(object$drmList[[1]]$fct$name == "LN.3u"){
#         constant <- qnorm(z)
#         closedExpression <- paste0("exp(",constant,"/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "LN.4"){
#         constant <- qnorm(bmr)
#         closedExpression <- paste0("exp(",constant,"/b)*e") # this is wrong
#       }
#       # weibull 1
#       else if(object$drmList[[1]]$fct$name == "W1.2"){
#         closedExpression <- paste0("(-log(",z,"))^(1/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "W1.3"){
#         closedExpression <- paste0("(-log(",z,"/d))^(1/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "W1.3u"){
#         closedExpression <- paste0("(-log((",z,"-c)/(1-c)))^(1/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "W1.4"){
#         closedExpression <- paste0("(-log((",z,"-c)/(d-c)))^(1/b)*e")
#       } 
#       # weibull 2
#       else if(object$drmList[[1]]$fct$name == "W2.2"){
#         closedExpression <- paste0("(-log(1-",z,"))^(1/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "W2.3"){
#         closedExpression <- paste0("(-log(1-",z,"/d))^(1/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "W2.3u"){
#         closedExpression <- paste0("(-log(1-(",z,"-c)/(1-c)))^(1/b)*e")
#       } else if(object$drmList[[1]]$fct$name == "W2.4"){
#         closedExpression <- paste0("(-log(1-(",z,"-c)/(d-c)))^(1/b)*e")
#       } 
#     }
#     BMD.pooled <- mjust(object$drmList,
#                         as.list(rep(closedExpression, length(object$drmList))),
#                         seType = "san")
#     tmp <- confint(multcomp:::glht(multcomp:::parm(BMD.pooled[["coef"]][,1],
#                                                    BMD.pooled[["covar"]]), linfct = matrix(rep(1/length(object$drmList), length(object$drmList)),1,length(object$drmList))),
#                    level = 1 - (1-level)*2)
#     CI <- tmp$confint[1,2:3]
#   }
#   CI
# }
