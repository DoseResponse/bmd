# drc version 3.2-0

test_that("bmd profile interval for ryegrass model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", level = 0.95, display = FALSE, interval = "profile")
  
  bmr <- 0.1
  def <- "extra"
  backgType <- "modelBased"
  interval <- "profile"
  backg <- NA
  controlSD <- NA
  respTrans <- "none"
  level <- 0.9

  
  # model unchanged
  expect_equal(as.numeric(logLik(object0)), -16.1551432181693)
  expect_equal(unname(coef(object0)), c(2.98221907127508,0.481413188446126,7.7929582936999,3.0579549665481))
  
  # profile interval calculation
  bmrScaled <- result$bmrScaled
  tmpVals <- ED(object0, bmrScaled, interval = "delta",
                level = level, type = "absolute", vcov. = vcov, display = FALSE)[,c("Lower", "Upper"), drop = FALSE]
  # expect_equal(unname(drop(tmpVals)), c(1.14157927967681,1.78583203136398)) # THIS MIGHT CHANGE BETWEEN DRC VERSIONS
  
  # tmpVals[,] <- c(1.14157927967681,1.78583203136398)
  # The mystery is resolved: The change in confidence intervals is due to different tmpVals. 
  # Tolerance in uniroot in bmdProfileCI has been lowered to 1e-8 to circumvent this issue in the future. 
  
  profileGridSize <- 20
  slope <- drop(ifelse(object0$curve[[1]](0)-object0$curve[[1]](Inf)>0,"decreasing","increasing"))
  expect_equal(slope, "decreasing")
  
  tmpVals[,"Lower"] <- max(tmpVals[,"Lower"], 1e-16) # ensure positive lower limit is supplied to bmdProfileCI
  
  bmdVal <- result$Results[1,1]
  tmpInterval <- bmdProfileCI(object0, slope, bmr, backgType, backg, def, respTrans, level = level, gridSize = profileGridSize,
                              bmdEst = bmdVal, lower = tmpVals[,"Lower"], upper = tmpVals[,"Upper"])
  expect_equal(unname(tmpInterval), c(1.19111858748113, 1.77626367708709)) # c(1.19113157178987,1.77626232094376)
  
  # go into bmdProfileCI function
  n <- object0$sumList$lenData
  dose <- object0$dataList$dose
  response <- object0$dataList$resp
  start <- coef(object0)[-length(coef(object0))]
  lower = tmpVals[,"Lower"]
  upper = tmpVals[,"Upper"]
  
  curveRepar <- getCurveRepar(object0, slope, bmr, backgType, backg, #controlSD, 
                              def, respTrans)
  profileLogLikFixedBmd <- getProfileLogLikFixedBmd(object0, curveRepar, bmr, start)
  expect_equal(profileLogLikFixedBmd(1.44), -16.1649093136363) # consistent logLikelihood for reparametrised function in arbitrary value
  expect_equal(profileLogLikFixedBmd(bmdVal), as.numeric(logLik(object0))) # consistent logLikelihood for reparametrised function and original function
  
  quant <- qchisq(p = level, df = 1)
  llMod <- as.numeric(logLik(object0))
  
  grid <- seq(lower, upper, length.out = profileGridSize)
  grid <- sort(c(grid, bmdVal)) # adding bmd estimate to ensure we have at least one grid point where H0 is accepted
  
  llVals <- sapply(grid, profileLogLikFixedBmd)
  accept <-  2 * (llMod - llVals) <= quant
  
  # Then, search for endpoints of CI between grid points
  CIlower <- uniroot(function(x) 2 * (llMod - profileLogLikFixedBmd(x)) - quant,
                     lower = grid[which(grid == min(grid[accept])) - 1],
                     upper = min(grid[accept]), tol = 1e-8)$root |> try(silent = TRUE) |> as.numeric()
  CIupper <- uniroot(function(x) 2 * (llMod - profileLogLikFixedBmd(x)) - quant,
                     lower = max(grid[accept]),
                     upper = grid[which(grid == max(grid[accept])) + 1], tol = 1e-8)$root |> try(silent = TRUE) |> as.numeric()
  
  expect_equal(CIlower, result$interval[1])
  expect_equal(CIlower, unname(tmpInterval[1]))
  expect_equal(CIlower, 1.19111858748113) # 1.19113157178987
  expect_equal(CIupper, result$interval[2])
  expect_equal(CIupper, unname(tmpInterval[2]))
  expect_equal(CIupper, 1.77626367708709) # 1.77626232094376
})
