.tukeytrendfit <- function (y, x, scaling = c("ari", "ord", "arilog"), ctype = NULL, ddf = c("residual", 
                                                                                                "KR", "PB"), d0shift = 1) 
{
  fit <- lm(y ~ x)
  # dose <- x
  ddf <- "residual"
  if(min(x) > 0){
    x.log <- x
    arilog <- log
  } else {
    x.log <- x - min(x)
    d0shift = d0shift
    arilog <- function(z){
      if(z == 0){
        x.unique <- sort(unique(x.log))
        log(x.unique[2]) - d0shift * (x.unique[2] - x.unique[1]) / (x.unique[3] - x.unique[2]) * (log(x.unique[3]) - log(x.unique[2]))
      } else
        log(z)
    }
    arilog <- Vectorize(arilog)
  }
  
  DAT <- data.frame(x=x,y=y)
  TDAT <- cbind(DAT,
                    xari=x,
                    xord=as.numeric(factor(x))-1,
                    xarilog=arilog(x.log))
  TNAM <- colnames(TDAT)[-(1:2)]
  SCAL <- scaling
  
  MLIST <- list()
  for (i in seq(along.with = SCAL)) {
    FORMI <- as.formula(paste(". ~ . - x + ", TNAM[i], 
                              sep = ""))
    MLIST[[i]] <- update(fit, FORMI, data = TDAT, na.action = "na.exclude")
  }
  names(MLIST) <- TNAM
  
  MMM <- MLIST
  class(MMM) <- "mmm"
  
  MLF <- as.list(paste(TNAM, " = 0", sep = ""))
  for (i in 1:3) {
    MLF[[i]] <- multcomp::glht(model = MMM[[i]], linfct = MLF[[i]])$linfct
  }
  names(MLF) <- TNAM
  class(MLF) <- "mlf"
  DF <- unlist(lapply(MLIST, df.residual))
  
  CALL <- fit$call
  cCALL <- as.character(CALL)
  cCALL[2] <- strsplit(cCALL[2], split = "[ ~ ]")[[1]][1]
  INFO <- paste(cCALL[1:2], collapse = ".")
  MODINFO <- list(modelinfo = INFO, initcall = CALL)
  
  # MODINFO <- getmodelinfo(fit)
  OUT <- c(list(mmm = MMM, mlf = MLF, df = DF), MODINFO)
  class(OUT) <- "tukeytrend"
  return(OUT)
}