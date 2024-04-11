bmdProfileCI <- function(object, slope, bmr, backgType, backg, # controlSD, 
                         def, respTrans, start, level = 0.95, gridSize = 10, bmdEst, lower, upper){
  n <- object$sumList$lenData
  dose <- object$dataList$dose
  response <- object$dataList$resp
  if(missing(start)){
    start <- coef(object)[-length(coef(object))]
  }
  
  curveRepar <- getCurveRepar(object, slope, bmr, backgType, backg, #controlSD, 
                              def, respTrans)
  profileLogLikFixedBmd <- getProfileLogLikFixedBmd(object, curveRepar, bmr, start)
  
  quant <- qchisq(p = level, df = 1)
  llMod <- as.numeric(logLik(object))
  
  # lower and upper bounds
  accept0 <- 2 * (llMod - profileLogLikFixedBmd(lower)) <= quant
  while(accept0){
    lower <- lower/2
    accept0 <-  2 * (llMod - profileLogLikFixedBmd(lower)) <= quant
  }
  
  accept0 <- 2 * (llMod - profileLogLikFixedBmd(upper)) <= quant
  while(accept0){
    upper <- upper*2
    accept0 <-  2 * (llMod - profileLogLikFixedBmd(upper)) <= quant
  }
  
  # Grid search
  grid <- seq(lower, upper, length.out = gridSize)
  grid <- sort(c(grid, bmdEst)) # adding bmd estimate to ensure we have at least one grid point where H0 is accepted
  llVals <- sapply(grid, profileLogLikFixedBmd)
  accept <-  2 * (llMod - llVals) <= quant
  
  # Then, search for endpoints of CI between grid points
  CIlower <- uniroot(function(x) 2 * (llMod - profileLogLikFixedBmd(x)) - quant,
                     lower = grid[which(grid == min(grid[accept])) - 1],
                     upper = min(grid[accept]))$root |> try(silent = TRUE) |> as.numeric()
  CIupper <- uniroot(function(x) 2 * (llMod - profileLogLikFixedBmd(x)) - quant,
                     lower = max(grid[accept]),
                     upper = grid[which(grid == max(grid[accept])) + 1])$root |> try(silent = TRUE) |> as.numeric()
  
  c(BMDL = CIlower, BMDU = CIupper)
}