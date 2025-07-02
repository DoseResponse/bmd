# Tests for bmdMA function
# - Arguments and structure
#   - Missing arguments
#   - Correct backgType and def accepted
#   - modelWeights argument
# - Simple model 
#   - correct bmd estimate for all (all definitions)
# - Ryegrass model (continuous)
#   - correct bmd estimate (all definitions)
# - TCDD model (binomial)
#   - correct bmd estimate (excess + additional)
# - S.alba model (continuous with multiple curves)
#   - correct bmd estimate (point, extra, hybridExc)
# - Decreasing binomial model with multiple curves
#   - correct bmd estimate (point, extra, hybridExc)

# Arguments and structure -------------------------------------------------

test_that("bmdMA function handles missing required arguments", {
  object0 <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2())
  modelList0 <- list(object0)
  
  expect_error(bmdMA(), 'argument "modelList" is missing, with no default')
  expect_error(bmdMA(modelList0), "def is missing")
  expect_error(bmdMA(list(lm(1:10 ~ 1))), 'modelList must be a list of models of class "drc"')
  expect_error(bmdMA(modelList0, def = "relative", backgType = "modelBased", bmr = 0.1), 'argument "modelWeights" is missing, with no default')
  expect_error(bmdMA(modelList0, modelWeights = "invalid_weights", def = "relative", backgType = "modelBased", bmr = 0.1), 'modelWeights must either be "AIC", "BIC", "Stack", "Stacking" or a numeric vector of same length as modelList')
  expect_error(bmdMA(modelList0, modelWeights = c(0.1,0.9), def = "relative", backgType = "modelBased", bmr = 0.1), 'modelWeights must either be "AIC", "BIC", "Stack", "Stacking" or a numeric vector of same length as modelList')
  expect_error(bmdMA(modelList0, modelWeights = "AIC", def = "relative", backgType = "modelBased", bmr = 0.1), 'Specify model averaging type. Options are "curve", "bootstrap", "Kang" and "Buckland"')
  expect_error(bmdMA(modelList0, modelWeights = "AIC", def = "invalid_def", backgType = "modelBased", type = "Kang"), "Could not recognize def")
  expect_error(bmdMA(modelList0, modelWeights = "AIC", def = "relative", backgType = "invalid_type", type = "Kang"), "Could not recognize backgType")
  expect_error(bmdMA(modelList0, modelWeights = "AIC", def = "relative", backgType = "modelBased", type = "Kang"), 'argument "bmr" is missing, with no default')
  expect_error(bmdMA(modelList0, modelWeights = "AIC", def = "relative", backgType = "modelBased", type = "Kang", bootstrapType = "invalid_type"), '"bootstrapType" not recognised. Options are: "nonparametric" and "parametric"')
})


test_that("bmdMA function accepts correct def", {
  object.cont <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2())
  object.binom <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2(), type = "binomial")
  object.poisson <- drm(y ~ x, data = data.frame(x = 1:5, y = c(12,11,3,0,0)), fct = LL.3(), type = "Poisson")
  
  # Binomial bmd def with continuous model
  expect_error(bmdMA(list(object.cont), modelWeights = 1, def = "excess", backgType = "modelBased", type = "Kang"), '"excess" is not available for continuous data')
  expect_error(bmdMA(list(object.cont), modelWeights = 1, def = "additional", backgType = "modelBased", type = "Kang"), '"additional" is not available for continuous data')
  
  # Binomial bmd def with Poisson model
  expect_error(bmdMA(list(object.poisson), modelWeights = 1, def = "excess", backgType = "modelBased", type = "Kang"), '"excess" is not available for count data')
  expect_error(bmdMA(list(object.poisson), modelWeights = 1, def = "additional", backgType = "modelBased", type = "Kang"), '"additional" is not available for count data')
  
  # Cont bmd def with binomial model
  expect_error(bmdMA(list(object.binom), modelWeights = 1, def = "relative", backgType = "modelBased", type = "Kang"), '"relative" is not available for quantal data')
  expect_error(bmdMA(list(object.binom), modelWeights = 1, def = "extra", backgType = "modelBased", type = "Kang"), '"extra" is not available for quantal data')
  expect_error(bmdMA(list(object.binom), modelWeights = 1, def = "added", backgType = "modelBased", type = "Kang"), '"added" is not available for quantal data')
  expect_error(bmdMA(list(object.binom), modelWeights = 1, def = "hybridExc", backgType = "modelBased", type = "Kang"), '"hybridExc" is not available for quantal data')
  expect_error(bmdMA(list(object.binom), modelWeights = 1, def = "hybridAdd", backgType = "modelBased", type = "Kang"), '"hybridAdd" is not available for quantal data')
})

test_that("bmdMA function returns expected structure", {
  object0 <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmdMA(list(object0), type = "Kang", modelWeights = "AIC", bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  
  expect_type(result, "list")
  expect_named(result, c("Results", "Boot.samples.used", "interval", "SE", "modelWeights"))
  expect_s3_class(result, "bmd")
})

test_that("bmdMA function handles modelWeights argument", {
  # data and fitted models
  data0 <- data.frame(x = rep(c(0,0.5, 1, 2, 4), 3), 
                      y = c(0.990994455038216, 0.743019106108585, 0.166411612682373, 0.0439546959424475, 0.183799321812153,
                            0.884930607478322, 1.02159286444501, 0.176423458342167, 0.108887606895092, 0.103881440788692, 
                            1.06023067243334, 0.921135238829789, 0.371858798338925, -0.0499408692067165, -0.204825758120983))
  object.LL.2 <- drm(y ~ x, data = data0, fct = LL.2()) 
  object.LN.2 <- drm(y ~ x, data = data0, fct = LN.2()) 
  object.W1.2 <- drm(y ~ x, data = data0, fct = W1.2()) 
  object.W2.2 <- drm(y ~ x, data = data0, fct = W2.2()) 
  modelList0 <- list(object.LL.2, object.LN.2, object.W1.2, object.W2.2)
  
  # bmd estimates on models
  bmdList0 <- lapply(modelList0, function(x) bmd(x, bmr = 0.08, backgType = "modelBased", def = "relative", display = FALSE))
  bmdVals <- sapply(bmdList0, function(x) x$Results[1,1])
  bmdlVals <- sapply(bmdList0, function(x) x$Results[1,2])
  
  # model weights
  manWeights0 <- c(0.3, 0.2, 0.2, 0.3)
  AICWeights0 <- exp(-1/2 * (sapply(modelList0, AIC) - min(sapply(modelList0, AIC)))) / sum(exp(-1/2 * (sapply(modelList0, AIC) - min(sapply(modelList0, AIC)))) )
  BICWeights0 <- exp(-1/2 * sapply(modelList0, BIC)) / sum(exp(-1/2 * sapply(modelList0, BIC)) )
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  stackingWeights0 <- getStackingWeights(modelList0, nSplits = 3)
  
  # apply bmdMA
  bmdMAManWeights <- bmdMA(modelList0, modelWeights = manWeights0, bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE)
  bmdMAAICWeights <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE)
  bmdMABICWeights <- bmdMA(modelList0, modelWeights = "BIC", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE)
  bmdMAStackingWeights <- bmdMA(modelList0, modelWeights = "Stack", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE, stackingSeed = 1, stackingSplits = 3)
  
  # checks
  expect_equal(bmdMAManWeights$modelWeights, manWeights0)
  expect_equal(bmdMAManWeights$Results[1,1], sum(bmdVals * manWeights0))
  expect_equal(bmdMAManWeights$Results[1,2], sum(bmdlVals * manWeights0))
  
  expect_equal(bmdMAAICWeights$modelWeights, AICWeights0)
  expect_equal(bmdMAAICWeights$Results[1,1], sum(bmdVals * AICWeights0))
  expect_equal(bmdMAAICWeights$Results[1,2], sum(bmdlVals * AICWeights0))
  
  expect_equal(bmdMABICWeights$modelWeights, BICWeights0)
  expect_equal(bmdMABICWeights$Results[1,1], sum(bmdVals * BICWeights0))
  expect_equal(bmdMABICWeights$Results[1,2], sum(bmdlVals * BICWeights0))
  
  expect_equal(bmdMAStackingWeights$modelWeights, stackingWeights0, tolerance = 1e-4)
  expect_equal(bmdMAStackingWeights$Results[1,1], sum(bmdVals * stackingWeights0), tolerance = 1e-4)
  expect_equal(bmdMAStackingWeights$Results[1,2], sum(bmdlVals * stackingWeights0), tolerance = 1e-4)
})


test_that("bmdMA function handles type argument", {
  # data and fitted models
  data0 <- data.frame(x = rep(c(0,0.5, 1, 2, 4), 3), 
                      y = c(0.990994455038216, 0.743019106108585, 0.166411612682373, 0.0439546959424475, 0.183799321812153,
                            0.884930607478322, 1.02159286444501, 0.176423458342167, 0.108887606895092, 0.103881440788692, 
                            1.06023067243334, 0.921135238829789, 0.371858798338925, -0.0499408692067165, -0.204825758120983))
  object.LL.2 <- drm(y ~ x, data = data0, fct = LL.2()) 
  object.LN.2 <- drm(y ~ x, data = data0, fct = LN.2()) 
  object.W1.2 <- drm(y ~ x, data = data0, fct = W1.2()) 
  object.W2.2 <- drm(y ~ x, data = data0, fct = W2.2()) 
  object.FPL.4 <- drm(y ~ x, data = data0, fct = FPL.4(p1 = -1, p2 = 1)) 
  modelList0 <- list(object.LL.2, object.LN.2, object.W1.2, object.W2.2, object.FPL.4)
  
  # bmd estimates on models
  bmdList0 <- lapply(modelList0, function(x) bmd(x, bmr = 0.08, backgType = "modelBased", def = "relative", display = FALSE))
  bmdVals <- sapply(bmdList0, function(x) x$Results[1,1])
  bmrScaledVals <- sapply(bmdList0, function(x) x$bmrScaled[1,1])
  bmdlVals <- sapply(bmdList0, function(x) x$interval[1,1])
  bmduVals <- sapply(bmdList0, function(x) x$interval[1,2])
  
  # model weights
  manWeights0 <- c(0.3, 0.18, 0.15, 0.22, 0.15)
  
  # apply bmdMA
  bmdMAKang <- bmdMA(modelList0, modelWeights = manWeights0, bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE)
  bmdMABuckland <- bmdMA(modelList0, modelWeights = manWeights0, bmr = 0.08, backgType = "modelBased", def = "relative", type = "Buckland", display = FALSE)
  set.seed(123)
  bmdMABoot <- suppressWarnings(bmdMA(modelList0, modelWeights = manWeights0, bmr = 0.08, backgType = "modelBased", def = "relative", type = "bootstrap", R = 100, display = FALSE, progressInfo = FALSE))
  set.seed(123)
  bmdMABootBCa <- suppressWarnings(bmdMA(modelList0, modelWeights = manWeights0, bmr = 0.08, backgType = "modelBased", def = "relative", type = "bootstrap", bootInterval = "BCa", R = 100, display = FALSE, progressInfo = FALSE))
  set.seed(123)
  bmdMACurve <- suppressWarnings(bmdMA(modelList0, modelWeights = manWeights0, bmr = 0.08, backgType = "modelBased", def = "relative", type = "curve", R = 100, display = FALSE, progressInfo = FALSE))
  set.seed(123)
  bmdMACurveBCa <- suppressWarnings(bmdMA(modelList0, modelWeights = manWeights0, bmr = 0.08, backgType = "modelBased", def = "relative", type = "curve", bootInterval = "BCa", R = 100, display = FALSE, progressInfo = FALSE))
  
  # checks
  expect_equal(bmdMAKang$Results[1,1], sum(manWeights0 * bmdVals))
  expect_equal(bmdMABuckland$Results[1,1], sum(manWeights0 * bmdVals))
  expect_equal(bmdMABoot$Results[1,1], sum(manWeights0 * bmdVals))
  expect_equal(bmdMABootBCa$Results[1,1], sum(manWeights0 * bmdVals))
  expect_equal(bmdMACurve$Results[1,1], 0.475485949896767) # manual calculation checked in v2.6.7
  expect_equal(bmdMACurveBCa$Results[1,1], 0.475485949896767) # manual calculation checked in v2.6.7
  expect_equal(sum(manWeights0 * bmrScaledVals), 
               sum(manWeights0 * sapply(modelList0, 
                                        function(x) x$curve[[1]](bmdMACurve$Results[1,1]))), 
               tolerance = 1e-4)
  
  # intervals
  expect_equal(c(bmdMAKang$interval[1,1], bmdMAKang$interval[1,2]),
               c(sum(manWeights0 * bmdlVals), sum(manWeights0 * bmduVals)), tolerance = 1e-6)
  expect_equal(c(bmdMABuckland$interval[1,1], bmdMABuckland$interval[1,2]),
               bmdMABuckland$Results[1,1] + qnorm(c(0.05, 0.95)) * bmdMABuckland$SE[1,1], tolerance = 1e-6)
  expect_equal(c(bmdMABoot$interval[1,1], bmdMABoot$interval[1,2]), 
               c(0.367163177724489, 0.63600316417734), tolerance = 1e-4) # manual calculation checked in v2.6.7
  expect_equal(bmdMABootBCa$Results[1,2], 0.386237633715302, tolerance = 1e-4)  # manual calculation checked in v2.6.7
  expect_equal(bmdMABootBCa$interval[1,2], "Not available for BCa bootstrap") 
  expect_equal(c(bmdMACurve$interval[1,1], bmdMACurve$interval[1,2]), 
               c(0.370341672214575, 0.622868602504527), tolerance = 1e-4) # manual calculation checked in v2.6.7
  expect_equal(bmdMACurveBCa$Results[1,2], 0.389400777586158, tolerance = 1e-4)  # manual calculation checked in v2.6.7
  expect_equal(bmdMABootBCa$interval[1,2], "Not available for BCa bootstrap") 
  
  # SE (only for Buckland)
  expect_equal(bmdMABuckland$SE[1, 1], 0.0785971536016299) # manual calculation checked in v2.6.7
})

test_that("bmdMA handles seed correctly when using Stacking weights",{
  # data and fitted models
  data0 <- data.frame(x = rep(c(0,0.5, 1, 2, 4), 3), 
                      y = c(0.990994455038216, 0.743019106108585, 0.166411612682373, 0.0439546959424475, 0.183799321812153,
                            0.884930607478322, 1.02159286444501, 0.176423458342167, 0.108887606895092, 0.103881440788692, 
                            1.06023067243334, 0.921135238829789, 0.371858798338925, -0.0499408692067165, -0.204825758120983))
  object.LL.2 <- drm(y ~ x, data = data0, fct = LL.2()) 
  object.LN.2 <- drm(y ~ x, data = data0, fct = LN.2()) 
  object.W1.2 <- drm(y ~ x, data = data0, fct = W1.2()) 
  object.W2.2 <- drm(y ~ x, data = data0, fct = W2.2()) 
  modelList0 <- list(object.LL.2, object.LN.2, object.W1.2, object.W2.2)
  
  # first seed of seed (1, 123)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  twoNormalVariables.seed1 <- rnorm(2)
  
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  oneNormalVariable.seed1 <- rnorm(1)
  bmdMAStackingWeights.seed123.inside <- bmdMA(modelList0, modelWeights = "Stack", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE, stackingSeed = 123, stackingSplits = 3)
  secondNormalVariable.seed1 <- rnorm(1)
  
  set.seed(123)
  bmdMAStackingWeights.seed123.outside <- bmdMA(modelList0, modelWeights = "Stack", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE, stackingSplits = 3)
  
  expect_equal(twoNormalVariables.seed1, c(oneNormalVariable.seed1, secondNormalVariable.seed1))
  expect_equal(bmdMAStackingWeights.seed123.inside$modelWeights, bmdMAStackingWeights.seed123.outside$modelWeights, tolerance = 1e-4)
  
  # second set of seed (156, 999)
  set.seed(156, kind = "Mersenne-Twister", normal.kind = "Inversion")
  twoNormalVariables.seed156 <- rnorm(2)
  
  set.seed(156, kind = "Mersenne-Twister", normal.kind = "Inversion")
  oneNormalVariable.seed156 <- rnorm(1)
  bmdMAStackingWeights.seed999.inside <- bmdMA(modelList0, modelWeights = "Stack", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE, stackingSeed = 999, stackingSplits = 3)
  secondNormalVariable.seed156 <- rnorm(1)
  
  set.seed(999, kind = "Mersenne-Twister", normal.kind = "Inversion")
  bmdMAStackingWeights.seed999.outside <- bmdMA(modelList0, modelWeights = "Stack", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE, stackingSplits = 3)
  
  expect_equal(twoNormalVariables.seed156, c(oneNormalVariable.seed156, secondNormalVariable.seed156))
  expect_equal(bmdMAStackingWeights.seed999.inside$modelWeights, bmdMAStackingWeights.seed999.outside$modelWeights, tolerance = 1e-4)
  
  # Check that the results are indeed different
  expect(all(twoNormalVariables.seed1 != twoNormalVariables.seed156), failure_message = "Variables from different seeds are not different")
  expect(all(c(oneNormalVariable.seed1, secondNormalVariable.seed1) != c(oneNormalVariable.seed156, secondNormalVariable.seed156)), 
         failure_message = "Variables from different seeds are not different")
  expect(sum(abs(bmdMAStackingWeights.seed123.inside$modelWeights - bmdMAStackingWeights.seed999.inside$modelWeights)) > 0.01, 
         failure_message = "Stacking Weights from different seeds are not different")
  expect(sum(abs(bmdMAStackingWeights.seed123.outside$modelWeights - bmdMAStackingWeights.seed999.outside$modelWeights)) > 0.01, 
         failure_message = "Stacking Weights from different seeds are not different")
})

test_that("bmdMA function output remains consistent", {
  # data and fitted models
  data0 <- data.frame(x = rep(c(0,0.5, 1, 2, 4), 3), 
                      y = c(0.990994455038216, 0.743019106108585, 0.166411612682373, 0.0439546959424475, 0.183799321812153,
                            0.884930607478322, 1.02159286444501, 0.176423458342167, 0.108887606895092, 0.103881440788692, 
                            1.06023067243334, 0.921135238829789, 0.371858798338925, -0.0499408692067165, -0.204825758120983))
  object.LL.2 <- drm(y ~ x, data = data0, fct = LL.2()) 
  object.LN.2 <- drm(y ~ x, data = data0, fct = LN.2()) 
  object.W1.2 <- drm(y ~ x, data = data0, fct = W1.2()) 
  object.W2.2 <- drm(y ~ x, data = data0, fct = W2.2()) 
  modelList0 <- list(object.LL.2, object.LN.2, object.W1.2, object.W2.2)
  
  bmdMAresult <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Buckland", display = FALSE)
  
  snapshot_data <- list(
    Results = as.list(bmdMAresult$Results),
    Boot.samples.used = as.list(bmdMAresult$Boot.samples.used),
    interval = as.list(bmdMAresult$interval),
    SE = as.list(bmdMAresult$SE)
  )
  
  # Store a snapshot of the entire result object
  expect_snapshot_value(snapshot_data, style = "deparse")
})




# ryegrass models ---------------------------------------------------------


test_that("bmdMA function computes BMD (point) correctly for ryegrass models", {
  # data and fitted models
  data0 <- drcData::ryegrass
  object.LL <- drm(rootl ~ conc, data = data0, fct = LL.4()) 
  object.LN <- drm(rootl ~ conc, data = data0, fct = LN.4()) 
  object.W1 <- drm(rootl ~ conc, data = data0, fct = W1.4()) 
  object.W2 <- drm(rootl ~ conc, data = data0, fct = W2.4()) 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultKang <- bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "Kang", display = FALSE)
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "Buckland", display = FALSE)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultBoot <- bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "bootstrap", R = 50, display = FALSE, progressInfo = FALSE)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultCurve <- bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "curve", R = 50, display = FALSE, progressInfo = FALSE)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultBootBCa <- bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "bootstrap", bootInterval = "BCa", R = 50, display = FALSE, progressInfo = FALSE)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultCurveBCa <- bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "curve", bootInterval = "BCa", R = 50, display = FALSE, progressInfo = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # Kang
  expect_true(!is.na(resultKang$Results[1, "BMD_MA"]))
  expect_equal(resultKang$Results[1, "BMD_MA"], 3.62285507866639)
  expect_equal(resultKang$interval[1,], c(BMDL_MA = 3.26863980433835, BMDU_MA = 3.97707035299443), tolerance = 1e-6)
  
  # Buckland
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 3.62285507866639)
  expect_equal(resultBuckland$SE[1,1], 0.209800492692605, tolerance = 1e-4)
  expect_equal(resultBuckland$interval[1,], c(BMDL_MA = 3.27776397732476, BMDU_MA = 3.96794618000803), tolerance = 1e-6)
  
  # Boot
  expect_true(!is.na(resultBoot$Results[1, "BMD_MA"]))
  expect_equal(resultBoot$Boot.samples.used, 50)
  expect_equal(resultBoot$Results[1, "BMD_MA"], 3.62285507866639)
  expect_equal(resultBoot$interval[1,], c(BMDL_MA = 3.15637429751311, BMDU_MA = 4.13206808564422), tolerance = 1e-4)
  
  # Curve
  expect_true(!is.na(resultCurve$Results[1, "BMD_MA"]))
  expect_equal(resultCurve$Boot.samples.used, 50)
  expect_equal(resultCurve$Results[1, "BMD_MA"], 3.62594076228422)
  expect_equal(resultCurve$interval[1,], c(BMDL_MA = 3.1364748814186, BMDU_MA = 4.13187952141919), tolerance = 1e-4)
  
  # BootBCa
  expect_true(all(!is.na(resultBootBCa$Results[, "BMD_MA"])))
  expect_equal(unname(resultBootBCa$Results[, "BMD_MA"]), c(3.62285507866639))
  expect_equal(resultBootBCa$Boot.samples.used, 50)
  expect_equal(unname(resultBootBCa$Results[,"BMDL_MA"]), c(3.28847336953601), tolerance = 1e-4)
  expect_equal(unname(resultBootBCa$interval[,"BMDU_MA"]), c("Not available for BCa bootstrap"))
  
  # CurveBCa
  expect_true(all(!is.na(resultCurveBCa$Results[, "BMD_MA"])))
  expect_equal(unname(resultCurveBCa$Results[, "BMD_MA"]), c(3.62594076228422))
  expect_equal(resultCurveBCa$Boot.samples.used, 50)
  expect_equal(unname(resultCurveBCa$Results[,"BMDL_MA"]), c(3.34556336964975), tolerance = 1e-4)
  expect_equal(unname(resultCurveBCa$interval[,"BMDU_MA"]), c("Not available for BCa bootstrap"))
})

test_that("bmdMA function computes BMD (extra) correctly for ryegrass models", {
  # data and fitted models
  data0 <- drcData::ryegrass
  object.LL <- drm(rootl ~ conc, data = data0, fct = LL.4()) 
  object.LN <- drm(rootl ~ conc, data = data0, fct = LN.4()) 
  object.W1 <- drm(rootl ~ conc, data = data0, fct = W1.4()) 
  object.W2 <- drm(rootl ~ conc, data = data0, fct = W2.4()) 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "extra", backgType = "modelBased", type = "Buckland", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 1.52703134707651)
  expect_equal(resultBuckland$SE[1,1], 0.210054516104069, tolerance = 1e-4)
  expect_equal(unname(resultBuckland$interval[1,]), c(1.1815224144052,1.87254027974782), tolerance = 1e-4)
})

test_that("bmdMA function computes BMD (relative) correctly for ryegrass models", {
  # data and fitted models
  data0 <- drcData::ryegrass
  object.LL <- drm(rootl ~ conc, data = data0, fct = LL.4()) 
  object.LN <- drm(rootl ~ conc, data = data0, fct = LN.4()) 
  object.W1 <- drm(rootl ~ conc, data = data0, fct = W1.4()) 
  object.W2 <- drm(rootl ~ conc, data = data0, fct = W2.4()) 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "relative", backgType = "modelBased", type = "Buckland", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 1.55622345982699)
  expect_equal(resultBuckland$SE[1,1], 0.189786051498867, tolerance = 1e-4)
  expect_equal(unname(resultBuckland$interval[1,]), c(1.24405318467428,1.8683937349797), tolerance = 1e-4)
})

test_that("bmdMA function computes BMD (added) correctly for ryegrass models", {
  # data and fitted models
  data0 <- drcData::ryegrass
  object.LL <- drm(rootl ~ conc, data = data0, fct = LL.4()) 
  object.LN <- drm(rootl ~ conc, data = data0, fct = LN.4()) 
  object.W1 <- drm(rootl ~ conc, data = data0, fct = W1.4()) 
  object.W2 <- drm(rootl ~ conc, data = data0, fct = W2.4()) 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "added", backgType = "modelBased", type = "Buckland", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 0.930418495860691)
  expect_equal(resultBuckland$SE[1,1], 0.264731500303042, tolerance = 1e-4)
  expect_equal(unname(resultBuckland$interval[1,]), c(0.494973927418928,1.36586306430245), tolerance = 1e-4)
})

test_that("bmdMA function computes BMD (hybridAdd with hybridSD background) correctly for ryegrass models", {
  # data and fitted models
  data0 <- drcData::ryegrass
  object.LL <- drm(rootl ~ conc, data = data0, fct = LL.4()) 
  object.LN <- drm(rootl ~ conc, data = data0, fct = LN.4()) 
  object.W1 <- drm(rootl ~ conc, data = data0, fct = W1.4()) 
  object.W2 <- drm(rootl ~ conc, data = data0, fct = W2.4()) 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, type = "Buckland", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 1.31942461589965)
  expect_equal(resultBuckland$SE[1,1], 0.215846392664226, tolerance = 1e-4)
  expect_equal(unname(resultBuckland$interval[1,]), c(0.964388894061504,1.67446033773779), tolerance = 1e-4)
})

test_that("bmdMA function computes BMD (hybridExc with hybridSD background) correctly for ryegrass models", {
  # data and fitted models
  data0 <- drcData::ryegrass
  object.LL <- drm(rootl ~ conc, data = data0, fct = LL.4()) 
  object.LN <- drm(rootl ~ conc, data = data0, fct = LN.4()) 
  object.W1 <- drm(rootl ~ conc, data = data0, fct = W1.4()) 
  object.W2 <- drm(rootl ~ conc, data = data0, fct = W2.4()) 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, type = "Buckland", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 1.31462740464376)
  expect_equal(resultBuckland$SE[1,1], 0.216356052118326, tolerance = 1e-4)
  expect_equal(unname(resultBuckland$interval[1,]), c(0.958753367604025,1.67050144168349), tolerance = 1e-4)
})

test_that("bmdMA function computes BMD (hybridExc with hybridPercentile background) correctly for ryegrass models", {
  # data and fitted models
  data0 <- drcData::ryegrass
  object.LL <- drm(rootl ~ conc, data = data0, fct = LL.4()) 
  object.LN <- drm(rootl ~ conc, data = data0, fct = LN.4()) 
  object.W1 <- drm(rootl ~ conc, data = data0, fct = W1.4()) 
  object.W2 <- drm(rootl ~ conc, data = data0, fct = W2.4()) 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, type = "Buckland", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 1.2022469808022)
  expect_equal(resultBuckland$SE[1,1], 0.229041579777732, tolerance = 1e-4)
  expect_equal(unname(resultBuckland$interval[1,]), c(0.825507107582098,1.57898685402229), tolerance = 1e-4)
})

test_that("bmdMA function computes BMD (relative) with log-transformed response correctly for ryegrass models", {
  # data and fitted models
  data0 <- drcData::ryegrass
  object.LL <- drm(log(rootl) ~ conc, data = data0, fct = LL.4()) 
  object.LN <- drm(log(rootl) ~ conc, data = data0, fct = LN.4()) 
  object.W1 <- drm(log(rootl) ~ conc, data = data0, fct = W1.4()) 
  object.W2 <- drm(log(rootl) ~ conc, data = data0, fct = W2.4()) 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", type = "Buckland", display = FALSE)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultBootstrap <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", type = "bootstrap", display = FALSE, progressInfo = FALSE)
  resultCurve <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", type = "curve", display = FALSE, progressInfo = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # resultBuckland
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 1.20542377217772)
  expect_equal(resultBuckland$SE[1,1], 0.353626557987876, tolerance = 1e-4)
  expect_equal(unname(resultBuckland$interval[1,]), c(0.623759845685001,1.78708769867044), tolerance = 1e-4)
  
  # resultBootstrap
  expect_true(all(!is.na(resultBootstrap$Results[, "BMD_MA"])))
  expect_equal(unname(resultBootstrap$Results[, "BMD_MA"]), c(1.20542377217772))
  expect_equal(resultBootstrap$Boot.samples.used, 1000)
  expect_equal(unname(resultBootstrap$interval[,"BMDL_MA"]), c(0.962987242762308), tolerance = 1e-4)
  expect_equal(unname(resultBootstrap$interval[,"BMDU_MA"]), c(1.66195885889725), tolerance = 1e-4)
  
  # resultCurve
  expect_true(all(!is.na(resultCurve$Results[, "BMD_MA"])))
  expect_equal(unname(resultCurve$Results[, "BMD_MA"]), c(1.2058581809568))
  expect_equal(resultCurve$Boot.samples.used, 1000)
  expect_equal(unname(resultCurve$interval[,"BMDL_MA"]), c(1.01029775850722), tolerance = 1e-4)
  expect_equal(unname(resultCurve$interval[,"BMDU_MA"]), c(1.75309602794422), tolerance = 1e-4)
})

test_that("bmdMA function computes BMD (relative) with square root-transformed response correctly for ryegrass models", {
  # data and fitted models
  data0 <- drcData::ryegrass
  object.LL.2 <- drm(sqrt(rootl) ~ conc, data = data0, fct = LL.4()) 
  object.LN.2 <- drm(sqrt(rootl) ~ conc, data = data0, fct = LN.4()) 
  object.W1.2 <- drm(sqrt(rootl) ~ conc, data = data0, fct = W1.4()) 
  object.W2.2 <- drm(sqrt(rootl) ~ conc, data = data0, fct = W2.4()) 
  modelList0 <- list(object.LL.2, object.LN.2, object.W1.2, object.W2.2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", type = "Buckland", display = FALSE)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultBootstrap <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", type = "bootstrap", display = FALSE, progressInfo = FALSE)
  resultCurve <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", type = "curve", display = FALSE, progressInfo = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # resultBuckland
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 1.50354728636659)
  expect_equal(resultBuckland$SE[1,1], 0.237926009341212, tolerance = 1e-4)
  expect_equal(unname(resultBuckland$interval[1,]), c(1.11219382695561,1.89490074577758), tolerance = 1e-4)
  
  # resultBootstrap
  expect_true(all(!is.na(resultBootstrap$Results[, "BMD_MA"])))
  expect_equal(unname(resultBootstrap$Results[, "BMD_MA"]), c(1.50354728636659))
  expect_equal(resultBootstrap$Boot.samples.used, 1000)
  expect_equal(unname(resultBootstrap$interval[,"BMDL_MA"]), c(1.34368562204519), tolerance = 1e-4)
  expect_equal(unname(resultBootstrap$interval[,"BMDU_MA"]), c(1.73002893650846), tolerance = 1e-4)
  
  # resultCurve
  expect_true(all(!is.na(resultCurve$Results[, "BMD_MA"])))
  expect_equal(unname(resultCurve$Results[, "BMD_MA"]), c(1.51713033393971))
  expect_equal(resultCurve$Boot.samples.used, 1000)
  expect_equal(unname(resultCurve$interval[,"BMDL_MA"]), c(1.36000313896168), tolerance = 1e-4)
  expect_equal(unname(resultCurve$interval[,"BMDU_MA"]), c(1.74646551140735), tolerance = 1e-4)
})






# TCDD models -----------------------------------------------------------


test_that("bmdMA function computes BMD (point) correctly for TCDD models", {
  # data and fitted models
  data0 <- drcData::TCDD
  object.LL <- drm(incidence/total ~ conc, weights = total, data = data0, fct = LL.4(), type = "binomial") 
  object.LN <- drm(incidence/total ~ conc, weights = total, data = data0, fct = LN.4(), type = "binomial") 
  object.W1 <- drm(incidence/total ~ conc, weights = total, data = data0, fct = W1.4(), type = "binomial") 
  object.W2 <- drm(incidence/total ~ conc, weights = total, data = data0, fct = W2.4(), type = "binomial") 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultKang <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.22, def = "point", backgType = "modelBased", type = "Kang", display = FALSE)
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.22, def = "point", backgType = "modelBased", type = "Buckland", display = FALSE)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultBoot <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 0.22, def = "point", backgType = "modelBased", type = "bootstrap", R = 50, display = FALSE, progressInfo = FALSE))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultBootBCa <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 0.22, def = "point", backgType = "modelBased", type = "bootstrap", bootInterval = "BCa", R = 50, display = FALSE, progressInfo = FALSE))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultCurve <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 0.22, def = "point", backgType = "modelBased", type = "curve", R = 50, display = FALSE, progressInfo = FALSE))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultCurveBCa <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 0.22, def = "point", backgType = "modelBased", type = "curve", bootInterval = "BCa", R = 50, display = FALSE, progressInfo = FALSE))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # Kang
  expect_true(!is.na(resultKang$Results[1, "BMD_MA"]))
  expect_equal(resultKang$Results[1, "BMD_MA"], 8.0976277220152)
  expect_equal(resultKang$interval[1,], c(BMDL_MA = -263.257027195816, BMDU_MA = 279.452282639846), tolerance = 1)
  
  # Buckland
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 8.0976277220152)
  expect_equal(resultBuckland$SE[1,1], 164.979529539839, tolerance = 1e-1)
  expect_equal(resultBuckland$interval[1,], c(BMDL_MA = 8.0976277220152 + qnorm(0.05) * 164.979529539839, BMDU_MA = 8.0976277220152 + qnorm(0.95) * 164.979529539839), tolerance = 1e-1)
  
  # Boot
  expect_true(!is.na(resultBoot$Results[1, "BMD_MA"]))
  expect_equal(resultBoot$Boot.samples.used, 48, tolerance = 1)
  expect_equal(resultBoot$Results[1, "BMD_MA"], 8.0976277220152)
  expect_equal(resultBoot$interval[1,], c(BMDL_MA = 6.27951876692377, BMDU_MA = 15.8918488967058), tolerance = 1e-1)
  
  # BootBCa
  expect_true(all(!is.na(resultBootBCa$Results[, "BMD_MA"])))
  expect_equal(unname(resultBootBCa$Results[, "BMD_MA"]), c(8.0976277220152))
  expect_equal(resultBootBCa$Boot.samples.used, 48, tolerance = 1)
  expect_equal(unname(resultBootBCa$Results[,"BMDL_MA"]), c(5.48916035257719), tolerance = 1e-1)
  expect_equal(unname(resultBootBCa$interval[,"BMDU_MA"]), c("Not available for BCa bootstrap"))
  
  # Curve
  expect_true(!is.na(resultCurve$Results[1, "BMD_MA"]))
  expect_equal(resultCurve$Boot.samples.used, 48, tolerance = 1)
  expect_equal(resultCurve$Results[1, "BMD_MA"], 7.93703689546043)
  expect_equal(resultCurve$interval[1,], c(BMDL_MA = 6.51501264650305, BMDU_MA = 16.4080707920735), tolerance = 1e-1)
  
  # CurveBCa
  expect_true(all(!is.na(resultCurveBCa$Results[, "BMD_MA"])))
  expect_equal(unname(resultCurveBCa$Results[, "BMD_MA"]), c(7.93703689546043))
  expect_equal(resultCurveBCa$Boot.samples.used, 48, tolerance = 1)
  expect_equal(unname(resultCurveBCa$Results[,"BMDL_MA"]), c(5.91280712553062), tolerance = 1e-1)
  expect_equal(unname(resultCurveBCa$interval[,"BMDU_MA"]), c("Not available for BCa bootstrap"))
})

test_that("bmdMA function computes BMD (excess) correctly for TCDD models", {
  # data and fitted models
  data0 <- drcData::TCDD
  object.LL <- drm(incidence/total ~ conc, weights = total, data = data0, fct = LL.4(), type = "binomial") 
  object.LN <- drm(incidence/total ~ conc, weights = total, data = data0, fct = LN.4(), type = "binomial") 
  object.W1 <- drm(incidence/total ~ conc, weights = total, data = data0, fct = W1.4(), type = "binomial") 
  object.W2 <- drm(incidence/total ~ conc, weights = total, data = data0, fct = W2.4(), type = "binomial") 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "excess", backgType = "modelBased", type = "Buckland", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 6.01418679961228)
  expect_equal(resultBuckland$SE[1,1], 141.305877641293, tolerance = 1e-1)
  expect_equal(unname(resultBuckland$interval[1,]), c(-226.41329854823,238.441672147454), tolerance = 1e-1)
})

test_that("bmdMA function computes BMD (additional) correctly for TCDD models", {
  # data and fitted models
  data0 <- drcData::TCDD
  object.LL <- drm(incidence/total ~ conc, weights = total, data = data0, fct = LL.4(), type = "binomial") 
  object.LN <- drm(incidence/total ~ conc, weights = total, data = data0, fct = LN.4(), type = "binomial") 
  object.W1 <- drm(incidence/total ~ conc, weights = total, data = data0, fct = W1.4(), type = "binomial") 
  object.W2 <- drm(incidence/total ~ conc, weights = total, data = data0, fct = W2.4(), type = "binomial") 
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "additional", backgType = "modelBased", type = "Buckland", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(resultBuckland$Results[1, "BMD_MA"]))
  expect_equal(resultBuckland$Results[1, "BMD_MA"], 6.05421250633682)
  expect_equal(resultBuckland$SE[1,1], 137.343568778351, tolerance = 1e-1)
  expect_equal(unname(resultBuckland$interval[1,]), c(-219.855854737192,231.964279749866), tolerance = 1e-1)
})




# S.alba models -----------------------------------------------------------

test_that("bmdMA function computes BMD (point) correctly for S.alba models", {
  # data and fitted models
  data0 <- drcData::S.alba
  object.LL <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LL.4())
  object.LN <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LN.4())
  object.W1 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W1.4())
  object.W2 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W2.4())
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultKang <- bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "Kang", display = FALSE)
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "Buckland", display = FALSE)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultBoot <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "bootstrap", R = 50, display = FALSE, progressInfo = FALSE))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultCurve <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "curve", R = 50, display = FALSE, progressInfo = FALSE))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultBootBCa <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "bootstrap", R = 50, bootInterval = "BCa", display = FALSE, progressInfo = FALSE))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultCurveBCa <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "curve", R = 50, bootInterval = "BCa", display = FALSE, progressInfo = FALSE))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # Kang
  expect_true(all(!is.na(resultKang$Results[, "BMD_MA"])))
  expect_equal(resultKang$Results[, "BMD_MA"], c(Glyphosate = 39.5922313568255, Bentazone = 22.1027508909314))
  expect_equal(resultKang$interval[,1], c(Glyphosate = 31.3181029781253, Bentazone = 18.8097460620025), tolerance = 1e-4)
  expect_equal(resultKang$interval[,2], c(Glyphosate = 47.8663597355256, Bentazone = 25.3957557198603), tolerance = 1e-4)
  
  # Buckland
  expect_true(all(!is.na(resultBuckland$Results[, "BMD_MA"])))
  expect_equal(resultBuckland$Results[, "BMD_MA"], c(Glyphosate = 39.5922313568255, Bentazone = 22.1027508909314))
  expect_equal(resultBuckland$SE[,1], c(Glyphosate = 4.96856818234909, Bentazone = 2.02139569706442), tolerance = 1e-4)
  expect_equal(resultBuckland$interval[,1], c(Glyphosate = 39.5922313568255, Bentazone = 22.1027508909314) +
                 qnorm(0.05) * c(Glyphosate = 4.96856818234909, Bentazone = 2.02139569706442), tolerance = 1e-4)
  expect_equal(resultBuckland$interval[,2], c(Glyphosate = 39.5922313568255, Bentazone = 22.1027508909314) +
                 qnorm(0.95) * c(Glyphosate = 4.96856818234909, Bentazone = 2.02139569706442), tolerance = 1e-4)
  
  # Boot
  expect_true(all(!is.na(resultBoot$Results[, "BMD_MA"])))
  expect_equal(resultBoot$Boot.samples.used, 50)
  expect_equal(resultBoot$Results[, "BMD_MA"], c(Glyphosate = 39.5922313568255, Bentazone = 22.1027508909314))
  expect_equal(resultBoot$interval[,1], c(Glyphosate = 29.4081958504551, Bentazone = 18.318530994358), tolerance = 1e-2)
  expect_equal(resultBoot$interval[,2], c(Glyphosate = 48.1571285773063, Bentazone = 26.2292827844959), tolerance = 1e-2)
  
  # Curve
  expect_true(all(!is.na(resultCurve$Results[, "BMD_MA"])))
  expect_equal(resultCurve$Boot.samples.used, 50)
  expect_equal(resultCurve$Results[, "BMD_MA"], c(Glyphosate = 39.6418457928274, Bentazone = 22.0090447085917))
  expect_equal(resultCurve$interval[,1], c(Glyphosate = 29.3925146204986, Bentazone = 18.283216864945), tolerance = 1e-2)
  expect_equal(resultCurve$interval[,2], c(Glyphosate = 47.0458800892314, Bentazone = 26.2673216228492), tolerance = 1e-2)
  
  # BootBCa
  expect_true(all(!is.na(resultBootBCa$Results[, "BMD_MA"])))
  expect_equal(unname(resultBootBCa$Results[, "BMD_MA"]), c(39.5922313568255, 22.1027508909314))
  expect_equal(resultBootBCa$Boot.samples.used, 50)
  expect_equal(unname(resultBootBCa$Results[,"BMDL_MA"]), c(31.0359096482404,18.5448190488073), tolerance = 1e-2)
  expect_equal(unname(resultBootBCa$interval[,"BMDU_MA"]), c("Not available for BCa bootstrap","Not available for BCa bootstrap"))
  
  # CurveBCa
  expect_true(all(!is.na(resultCurveBCa$Results[, "BMD_MA"])))
  expect_equal(unname(resultCurveBCa$Results[, "BMD_MA"]), c(39.6418457928274, 22.0090447085917))
  expect_equal(resultCurveBCa$Boot.samples.used, 50)
  expect_equal(unname(resultCurveBCa$Results[,"BMDL_MA"]), c(30.9098852016046,18.4608747446523), tolerance = 1e-2)
  expect_equal(unname(resultCurveBCa$interval[,"BMDU_MA"]), c("Not available for BCa bootstrap","Not available for BCa bootstrap"))
})

test_that("bmdMA function computes BMD (relative) correctly for S.alba models", {
  # data and fitted models
  data0 <- drcData::S.alba
  object.LL <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LL.4())
  object.LN <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LN.4())
  object.W1 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W1.4())
  object.W2 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W2.4())
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.08, def = "relative", backgType = "modelBased", type = "Buckland", display = FALSE)
  
  # Buckland
  expect_true(all(!is.na(resultBuckland$Results[, "BMD_MA"])))
  expect_equal(resultBuckland$Results[, "BMD_MA"], c(Glyphosate = 28.8918829608504, Bentazone = 19.0370848122154))
  expect_equal(resultBuckland$SE[,1], c(Glyphosate = 6.97566436363773, Bentazone = 2.43528735772882), tolerance = 1e-4)
  expect_equal(resultBuckland$interval[,1], c(Glyphosate = 28.8918829608504, Bentazone = 19.0370848122154) +
                 qnorm(0.05) * c(Glyphosate = 6.97566436363773, Bentazone = 2.43528735772882), tolerance = 1e-4)
  expect_equal(resultBuckland$interval[,2], c(Glyphosate = 28.8918829608504, Bentazone = 19.0370848122154) +
                 qnorm(0.95) * c(Glyphosate = 6.97566436363773, Bentazone = 2.43528735772882), tolerance = 1e-4)
})

test_that("bmdMA function computes BMD (hybridExc with hybridSD background) correctly for S.alba models", {
  # data and fitted models
  data0 <- drcData::S.alba
  object.LL <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LL.4())
  object.LN <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LN.4())
  object.W1 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W1.4())
  object.W2 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W2.4())
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, type = "Buckland", display = FALSE)
  
  # Buckland
  expect_true(all(!is.na(resultBuckland$Results[, "BMD_MA"])))
  expect_equal(resultBuckland$Results[, "BMD_MA"], c(Glyphosate = 28.878538008635, Bentazone = 19.097790315674))
  expect_equal(resultBuckland$SE[,1], c(Glyphosate = 7.15914360496773, Bentazone = 2.48849148981338), tolerance = 1e-4)
  expect_equal(resultBuckland$interval[,1], c(Glyphosate = 28.878538008635, Bentazone = 19.097790315674) +
                 qnorm(0.05) * c(Glyphosate = 7.15914360496773, Bentazone = 2.48849148981338), tolerance = 1e-4)
  expect_equal(resultBuckland$interval[,2], c(Glyphosate = 28.878538008635, Bentazone = 19.097790315674) +
                 qnorm(0.95) * c(Glyphosate = 7.15914360496773, Bentazone = 2.48849148981338), tolerance = 1e-4)
})

test_that("bmdMA function handles modelWeights argument on S.alba data with multiple curves", {
  # data and fitted models
  data0 <- drcData::S.alba
  object.LL <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LL.4())
  object.LN <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LN.4())
  object.W1 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W1.4())
  object.W2 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W2.4())
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # bmd estimates on models
  bmdList0 <- lapply(modelList0, function(x) bmd(x, bmr = 0.08, backgType = "modelBased", def = "relative", display = FALSE))
  bmdVals <- sapply(bmdList0, function(x) x$Results[,1])
  bmdlVals <- sapply(bmdList0, function(x) x$Results[,2])
  bmduVals <- sapply(bmdList0, function(x) x$interval[,2])
  
  # model weights
  manWeights0 <- c(0.3, 0.2, 0.2, 0.3)
  AICWeights0 <- exp(-1/2 * (sapply(modelList0, AIC) - min(sapply(modelList0, AIC)))) / sum(exp(-1/2 * (sapply(modelList0, AIC) - min(sapply(modelList0, AIC)))) )
  BICWeights0 <- exp(-1/2 * sapply(modelList0, BIC)) / sum(exp(-1/2 * sapply(modelList0, BIC)) )
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  stackingWeights0 <- getStackingWeights(modelList0, nSplits = 4)
  
  # apply bmdMA
  bmdMAManWeights <- bmdMA(modelList0, modelWeights = manWeights0, bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE)
  bmdMAAICWeights <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE)
  bmdMABICWeights <- bmdMA(modelList0, modelWeights = "BIC", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE)
  bmdMAStackingWeights <- bmdMA(modelList0, modelWeights = "Stack", bmr = 0.08, backgType = "modelBased", def = "relative", type = "Kang", display = FALSE, stackingSeed = 1, stackingSplits = 4)
  
  # checks
  expect_equal(bmdMAManWeights$modelWeights, manWeights0)
  expect_equal(bmdMAManWeights$Results[1,1], sum(bmdVals[1,] * manWeights0))
  expect_equal(bmdMAManWeights$Results[2,1], sum(bmdVals[2,] * manWeights0))
  expect_equal(bmdMAManWeights$Results[1,2], sum(bmdlVals[1,] * manWeights0))
  expect_equal(bmdMAManWeights$Results[2,2], sum(bmdlVals[2,] * manWeights0))
  expect_equal(bmdMAManWeights$interval[1,2], sum(bmduVals[1,] * manWeights0))
  expect_equal(bmdMAManWeights$interval[2,2], sum(bmduVals[2,] * manWeights0))
  
  expect_equal(bmdMAAICWeights$modelWeights, AICWeights0)
  expect_equal(bmdMAAICWeights$Results[1,1], sum(bmdVals[1,] * AICWeights0))
  expect_equal(bmdMAAICWeights$Results[2,1], sum(bmdVals[2,] * AICWeights0))
  expect_equal(bmdMAAICWeights$Results[1,2], sum(bmdlVals[1,] * AICWeights0))
  expect_equal(bmdMAAICWeights$Results[2,2], sum(bmdlVals[2,] * AICWeights0))
  expect_equal(bmdMAAICWeights$interval[1,2], sum(bmduVals[1,] * AICWeights0))
  expect_equal(bmdMAAICWeights$interval[2,2], sum(bmduVals[2,] * AICWeights0))
  
  expect_equal(bmdMABICWeights$modelWeights, BICWeights0)
  expect_equal(bmdMABICWeights$Results[1,1], sum(bmdVals[1,] * BICWeights0))
  expect_equal(bmdMABICWeights$Results[2,1], sum(bmdVals[2,] * BICWeights0))
  expect_equal(bmdMABICWeights$Results[1,2], sum(bmdlVals[1,] * BICWeights0))
  expect_equal(bmdMABICWeights$Results[2,2], sum(bmdlVals[2,] * BICWeights0))
  expect_equal(bmdMABICWeights$interval[1,2], sum(bmduVals[1,] * BICWeights0))
  expect_equal(bmdMABICWeights$interval[2,2], sum(bmduVals[2,] * BICWeights0))
  
  expect_equal(bmdMAStackingWeights$modelWeights, stackingWeights0)
  expect_equal(bmdMAStackingWeights$Results[1,1], sum(bmdVals[1,] * stackingWeights0), tolerance = 1e-4)
  expect_equal(bmdMAStackingWeights$Results[2,1], sum(bmdVals[2,] * stackingWeights0), tolerance = 1e-4)
  expect_equal(bmdMAStackingWeights$Results[1,2], sum(bmdlVals[1,] * stackingWeights0), tolerance = 1e-4)
  expect_equal(bmdMAStackingWeights$Results[2,2], sum(bmdlVals[2,] * stackingWeights0), tolerance = 1e-4)
  expect_equal(bmdMAStackingWeights$interval[1,2], sum(bmduVals[1,] * stackingWeights0), tolerance = 1e-4)
  expect_equal(bmdMAStackingWeights$interval[2,2], sum(bmduVals[2,] * stackingWeights0), tolerance = 1e-4)
})

test_that("bmdMA function computes BMD (relative) correctly for S.alba models fitted separately", {
  # data and fitted models
  data0 <- drcData::S.alba
  object.LL <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LL.4(), separate = TRUE)
  object.LN <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LN.4(), separate = TRUE)
  object.W1 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W1.4(), separate = TRUE)
  object.W2 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W2.4(), separate = TRUE)
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.08, def = "relative", backgType = "modelBased", type = "Buckland", display = FALSE)
  
  # Buckland
  expect_true(all(!is.na(resultBuckland$Results[, "BMD_MA"])))
  expect_equal(resultBuckland$Results[, "BMD_MA"], c(Glyphosate = 28.8997770186961, Bentazone = 19.0337863593213))
  expect_equal(resultBuckland$SE[,1], c(Glyphosate = 7.46204086748364, Bentazone = 2.24787601134509), tolerance = 1e-4)
  expect_equal(resultBuckland$interval[,1], c(Glyphosate = 28.8997770186961, Bentazone = 19.0337863593213) +
                 qnorm(0.05) * c(Glyphosate = 7.46204086748364, Bentazone = 2.24787601134509), tolerance = 1e-4)
  expect_equal(resultBuckland$interval[,2], c(Glyphosate = 28.8997770186961, Bentazone = 19.0337863593213) +
                 qnorm(0.95) * c(Glyphosate = 7.46204086748364, Bentazone = 2.24787601134509), tolerance = 1e-4)
  expect_equal(resultBuckland$modelWeights, matrix(c(0.289011738824658, 0.252712376289855, 
                                                     0.242970006217661, 0.238530168705929,
                                                     0.232912736014905, 0.254486486073599, 
                                                     0.235105518942777, 0.254270968930617), 
                                                   nrow = 2,
                                                   dimnames = list(c("Glyphosate", "Bentazone"), NULL)), tolerance = 1e-4)
})


test_that("bmdMA function output remains consistent with model with multiple curves", {
  # data and fitted models
  data0 <- drcData::S.alba
  object.LL <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LL.4())
  object.LN <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = LN.4())
  object.W1 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W1.4())
  object.W2 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = data0, fct = W2.4())
  modelList0 <- list(object.LL, object.LN, object.W1, object.W2)
  
  # results
  resultBuckland <- bmdMA(modelList0, modelWeights = "AIC", bmr = 0.08, def = "extra", backgType = "modelBased", type = "Buckland", display = FALSE)
  
  snapshot_data <- list(
    Results = as.list(resultBuckland$Results),
    Boot.samples.used = as.list(resultBuckland$Boot.samples.used),
    interval = as.list(resultBuckland$interval),
    SE = as.list(resultBuckland$SE)
  )
  
  # Store a snapshot of the entire result object
  expect_snapshot_value(snapshot_data, style = "deparse")
})






































# Decreasing binomial model with multiple curves --------------------------
test_that("bmdMA function computes BMD (point) correctly for Decreasing binomial model with multiple curves", {
  data0 <- data.frame(
    conc = c(0, rep(c(20, 40, 80, 160, 320), 2)),
    total = c(50, rep(20, 5*2)),
    alive = c(47, 14, 11, 6, 9, 6, 19, 11, 8, 5, 3),
    treat = c("Control", rep(c("Treat1", "Treat2"), each = 5))
  )
  
  object0.LL.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = LL.4(), type = "binomial")
  object0.LN.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = LN.4(), type = "binomial")
  object0.W1.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = W1.4(), type = "binomial")
  object0.W2.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = W2.4(), type = "binomial")
  modelList0 <- list(object0.LL.4, object0.LN.4, object0.W1.4, object0.W2.4)
  
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  invisible(capture.output({
    resultKang <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 0.77, def = "point", backgType = "modelBased", type = "Kang", display = FALSE))
    resultBuckland <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 0.77, def = "point", backgType = "modelBased", type = "Buckland", display = FALSE))
    resultBoot <- suppressWarnings(bmdMA(modelList0, modelWeights = "BIC", bmr = 0.77, def = "point", backgType = "modelBased", type = "bootstrap", R = 50, progressInfo = FALSE, display = FALSE))
    resultCurve <- suppressWarnings(bmdMA(modelList0, modelWeights = "BIC", bmr = 0.77, def = "point", backgType = "modelBased", type = "curve", R = 50, progressInfo = FALSE, display = FALSE))
    resultBootBCa <- suppressWarnings(bmdMA(modelList0, modelWeights = "BIC", bmr = 0.77, def = "point", backgType = "modelBased", type = "bootstrap", bootInterval = "BCa", R = 50, progressInfo = FALSE, display = FALSE))
    resultCurveBCa <- suppressWarnings(bmdMA(modelList0, modelWeights = "BIC", bmr = 0.77, def = "point", backgType = "modelBased", type = "curve", bootInterval = "BCa", R = 50, progressInfo = FALSE, display = FALSE))
  }))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # resultKang
  expect_true(all(!is.na(resultKang$Results[, "BMD_MA"])))
  expect_equal(unname(resultKang$Results[, "BMD_MA"]), c(16.4767601484833, 30.3381343833517))
  expect_equal(resultKang$Boot.samples.used, NA)
  expect_equal(unname(resultKang$interval[,"BMDL_MA"]), c(5.4148182829762,21.1390036616694), tolerance = 1)
  expect_equal(unname(resultKang$interval[,"BMDU_MA"]), c(27.5387020139904,39.5372651050341), tolerance = 1)
  # resultBuckland
  expect_true(all(!is.na(resultBuckland$Results[, "BMD_MA"])))
  expect_equal(unname(resultBuckland$Results[, "BMD_MA"]), c(16.4767601484833, 30.3381343833517))
  expect_equal(resultBuckland$Boot.samples.used, NA)
  expect_equal(unname(resultBuckland$interval[,"BMDL_MA"]), c(5.40824166637437,21.1373086792919), tolerance = 1e-1)
  expect_equal(unname(resultBuckland$interval[,"BMDU_MA"]), c(27.5452786305922,39.5389600874116), tolerance = 1e-1)
  # resultBoot
  expect_true(all(!is.na(resultBoot$Results[, "BMD_MA"])))
  expect_equal(unname(resultBoot$Results[, "BMD_MA"]), c(16.4767601484833, 30.3381343833517))
  expect_equal(resultBoot$Boot.samples.used, 50, tolerance = 1)
  expect_equal(unname(resultBoot$interval[,"BMDL_MA"]), c(1.55114527115658,24.7976584352616), tolerance = 1e-1)
  expect_equal(unname(resultBoot$interval[,"BMDU_MA"]), c(28.447104103199,39.2750258645232), tolerance = 1e-1)
  # resultCurve
  expect_true(all(!is.na(resultCurve$Results[, "BMD_MA"])))
  expect_equal(unname(resultCurve$Results[, "BMD_MA"]), c(16.5073556830358, 30.3309353075015))
  expect_equal(resultCurve$Boot.samples.used, 50, tolerance = 1)
  expect_equal(unname(resultCurve$interval[,"BMDL_MA"]), c(1.45044325358349,26.3872029783547), tolerance = 1e-1)
  expect_equal(unname(resultCurve$interval[,"BMDU_MA"]), c(23.4835755264987,39.9141782740096), tolerance = 1e-1)
  # resultBootBCa
  expect_true(all(!is.na(resultBootBCa$Results[, "BMD_MA"])))
  expect_equal(unname(resultBootBCa$Results[, "BMD_MA"]), c(16.4767601484833, 30.3381343833517))
  expect_equal(resultBootBCa$Boot.samples.used, 49, tolerance = 1)
  expect_equal(unname(resultBootBCa$Results[,"BMDL_MA"]), c(3.3270653573855,16.8290145235522), tolerance = 1e-1)
  expect_equal(unname(resultBootBCa$interval[,"BMDU_MA"]), c("Not available for BCa bootstrap","Not available for BCa bootstrap"))
  # resultCurveBCa
  expect_true(all(!is.na(resultCurveBCa$Results[, "BMD_MA"])))
  expect_equal(unname(resultCurveBCa$Results[, "BMD_MA"]), c(16.5073556830358, 30.3309353075015))
  expect_equal(resultCurveBCa$Boot.samples.used, 50, tolerance = 1)
  expect_equal(unname(resultCurveBCa$Results[,"BMDL_MA"]), c(0.361186441608513,21.6147459158837), tolerance = 1e-1)
  expect_equal(unname(resultCurveBCa$interval[,"BMDU_MA"]), c("Not available for BCa bootstrap","Not available for BCa bootstrap"))
})

test_that("bmdMA function computes BMD (excess) correctly for Decreasing binomial model with multiple curves", {
  data0 <- data.frame(
    conc = c(0, rep(c(20, 40, 80, 160, 320), 2)),
    total = c(50, rep(20, 5*2)),
    alive = c(47, 14, 11, 6, 9, 6, 19, 11, 8, 5, 3),
    treat = c("Control", rep(c("Treat1", "Treat2"), each = 5))
  )
  
  object0.LL.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = LL.4(), type = "binomial")
  object0.LN.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = LN.4(), type = "binomial")
  object0.W1.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = W1.4(), type = "binomial")
  object0.W2.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = W2.4(), type = "binomial")
  modelList0 <- list(object0.LL.4, object0.LN.4, object0.W1.4, object0.W2.4)
  
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  invisible(capture.output({
    resultKang <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "excess", backgType = "modelBased", type = "Kang", display = FALSE))
    resultBuckland <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "excess", backgType = "modelBased", type = "Buckland", display = FALSE))
    resultBoot <- suppressWarnings(bmdMA(modelList0, modelWeights = "BIC", bmr = 0.1, def = "excess", backgType = "modelBased", type = "bootstrap", R = 50, progressInfo = FALSE, display = FALSE))
    resultCurve <- suppressWarnings(bmdMA(modelList0, modelWeights = "BIC", bmr = 0.1, def = "excess", backgType = "modelBased", type = "curve", R = 50, progressInfo = FALSE, display = FALSE))
    resultBootBCa <- suppressWarnings(bmdMA(modelList0, modelWeights = "BIC", bmr = 0.1, def = "excess", backgType = "modelBased", type = "bootstrap", R = 50, bootInterval = "BCa", progressInfo = FALSE, display = FALSE))
    resultCurveBCa <- suppressWarnings(bmdMA(modelList0, modelWeights = "BIC", bmr = 0.1, def = "excess", backgType = "modelBased", type = "curve", R = 50, bootInterval = "BCa", progressInfo = FALSE, display = FALSE))
  }))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # resultKang
  expect_true(all(!is.na(resultKang$Results[, "BMD_MA"])))
  expect_equal(unname(resultKang$Results[, "BMD_MA"]), c(11.7349943036545, 23.7533160131826))
  expect_equal(resultKang$Boot.samples.used, NA)
  expect_equal(unname(resultKang$interval[,"BMDL_MA"]), c(0.0079507187171709,13.7125460518548), tolerance = 1e-2)
  expect_equal(unname(resultKang$interval[,"BMDU_MA"]), c(23.4620378885919,33.7940859745104), tolerance = 1e-2)
  # resultBuckland
  expect_true(all(!is.na(resultBuckland$Results[, "BMD_MA"])))
  expect_equal(unname(resultBuckland$Results[, "BMD_MA"]), c(11.7349943036545, 23.7533160131826))
  expect_equal(resultBuckland$Boot.samples.used, NA)
  expect_equal(unname(resultBuckland$interval[,"BMDL_MA"]), c(-0.114127421954805,13.4685618858526), tolerance = 1e-2)
  expect_equal(unname(resultBuckland$interval[,"BMDU_MA"]), c(23.5841160292639,34.0380701405126), tolerance = 1e-2)
  # resultBoot
  expect_true(all(!is.na(resultBoot$Results[, "BMD_MA"])))
  expect_equal(unname(resultBoot$Results[, "BMD_MA"]), c(11.7349943036545, 23.7533160131826))
  expect_equal(resultBoot$Boot.samples.used, 50, tolerance = 1)
  expect_equal(unname(resultBoot$interval[,"BMDL_MA"]), c(0.176508789719701,7.73106039713397), tolerance = 1e-1)
  expect_equal(unname(resultBoot$interval[,"BMDU_MA"]), c(24.7982532099151,32.7273607506471), tolerance = 1e-1)
  # resultCurve
  expect_true(all(!is.na(resultCurve$Results[, "BMD_MA"])))
  expect_equal(unname(resultCurve$Results[, "BMD_MA"]), c(11.9000482326248, 24.0167032503988))
  expect_equal(resultCurve$Boot.samples.used, 49, tolerance = 1)
  expect_equal(unname(resultCurve$interval[,"BMDL_MA"]), c(0.212082383992232,17.245459902478), tolerance = 1e-1)
  expect_equal(unname(resultCurve$interval[,"BMDU_MA"]), c(19.1178603752347, 34.4447842377776), tolerance = 1e-1)
  # resultBootBCa
  expect_true(all(!is.na(resultBootBCa$Results[, "BMD_MA"])))
  expect_equal(unname(resultBootBCa$Results[, "BMD_MA"]), c(11.7349943036545, 23.7533160131826))
  expect_equal(resultBootBCa$Boot.samples.used, 49, tolerance = 1)
  expect_equal(unname(resultBootBCa$Results[,"BMDL_MA"]), c(0.422313711087127,2.45309962249564), tolerance = 1e-1)
  expect_equal(unname(resultBootBCa$interval[,"BMDU_MA"]), c("Not available for BCa bootstrap","Not available for BCa bootstrap"))
  # ResultCurveBCa
  expect_true(all(!is.na(resultCurveBCa$Results[, "BMD_MA"])))
  expect_equal(unname(resultCurveBCa$Results[, "BMD_MA"]), c(11.9000482326248, 24.0167032503988))
  expect_equal(resultCurveBCa$Boot.samples.used, 49, tolerance = 1)
  expect_equal(unname(resultCurveBCa$Results[,"BMDL_MA"]), c(1.01766053328384,11.5161237490573), tolerance = 1e-1)
  expect_equal(unname(resultCurveBCa$interval[,"BMDU_MA"]), c("Not available for BCa bootstrap","Not available for BCa bootstrap"))
})

test_that("bmdMA function computes BMD (additional) correctly for Decreasing binomial model with multiple curves", {
  data0 <- data.frame(
    conc = c(0, rep(c(20, 40, 80, 160, 320), 2)),
    total = c(50, rep(20, 5*2)),
    alive = c(47, 14, 11, 6, 9, 6, 19, 11, 8, 5, 3),
    treat = c("Control", rep(c("Treat1", "Treat2"), each = 5))
  )
  
  object0.LL.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = LL.4(), type = "binomial")
  object0.LN.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = LN.4(), type = "binomial")
  object0.W1.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = W1.4(), type = "binomial")
  object0.W2.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = W2.4(), type = "binomial")
  modelList0 <- list(object0.LL.4, object0.LN.4, object0.W1.4, object0.W2.4)
  
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  invisible(capture.output({
    resultKang <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "additional", backgType = "modelBased", type = "Kang", display = FALSE))
    resultBuckland <- suppressWarnings(bmdMA(modelList0, modelWeights = "AIC", bmr = 0.1, def = "additional", backgType = "modelBased", type = "Buckland", display = FALSE))
    resultBoot <- suppressWarnings(bmdMA(modelList0, modelWeights = "BIC", bmr = 0.1, def = "additional", backgType = "modelBased", type = "bootstrap", R = 50, progressInfo = FALSE, display = FALSE))
    resultCurve <- suppressWarnings(bmdMA(modelList0, modelWeights = "BIC", bmr = 0.1, def = "additional", backgType = "modelBased", type = "curve", R = 50, progressInfo = FALSE, display = FALSE))
  }))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # resultKang
  expect_true(all(!is.na(resultKang$Results[, "BMD_MA"])))
  expect_equal(unname(resultKang$Results[, "BMD_MA"]), c(12.0475135843298, 24.2168984689448))
  expect_equal(resultKang$Boot.samples.used, NA)
  expect_equal(unname(resultKang$interval[,"BMDL_MA"]), c(0.273646037712913,14.129549861443), tolerance = 1e-2)
  expect_equal(unname(resultKang$interval[,"BMDU_MA"]), c(23.8213811309467,34.3042470764466), tolerance = 1e-2)
  # resultBuckland
  expect_true(all(!is.na(resultBuckland$Results[, "BMD_MA"])))
  expect_equal(unname(resultBuckland$Results[, "BMD_MA"]), c(12.0475135843298, 24.2168984689448))
  expect_equal(resultBuckland$Boot.samples.used, NA)
  expect_equal(unname(resultBuckland$interval[,"BMDL_MA"]), c(0.164984994060264,13.9222959869599), tolerance = 1e-2)
  expect_equal(unname(resultBuckland$interval[,"BMDU_MA"]), c(23.9300421745994,34.5115009509297), tolerance = 1e-2)
  # resultBoot
  expect_true(all(!is.na(resultBoot$Results[, "BMD_MA"])))
  expect_equal(unname(resultBoot$Results[, "BMD_MA"]), c(12.0475135843298, 24.2168984689448))
  expect_equal(resultBoot$Boot.samples.used, 50, tolerance = 1)
  expect_equal(unname(resultBoot$interval[,"BMDL_MA"]), c(0.160028906689239,7.38437587559939), tolerance = 1e-1)
  expect_equal(unname(resultBoot$interval[,"BMDU_MA"]), c(25.2916075533545,33.2573011197188), tolerance = 1e-1)
  # resultCurve
  expect_true(all(!is.na(resultCurve$Results[, "BMD_MA"])))
  expect_equal(unname(resultCurve$Results[, "BMD_MA"]), c(12.2047330511244, 24.4616008187769))
  expect_equal(resultCurve$Boot.samples.used, 48, tolerance = 1)
  expect_equal(unname(resultCurve$interval[,"BMDL_MA"]), c(0.641340012245435,18.4858704146623), tolerance = 1e-1)
  expect_equal(unname(resultCurve$interval[,"BMDU_MA"]), c(19.5768964306862,34.7795760773764), tolerance = 1e-1)
})


test_that("bmdMA function computes BMD (point with stacking weights) correctly for decreasing binomial model with multiple curves", {
  data0 <- data.frame(
    conc = c(0, rep(c(20, 40, 80, 160, 320), 2)),
    total = c(50, rep(20, 5*2)),
    alive = c(47, 14, 11, 6, 9, 6, 19, 11, 8, 5, 3),
    treat = c("Control", rep(c("Treat1", "Treat2"), each = 5))
  )
  
  object0.LL.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = LL.4(), type = "binomial")
  object0.LN.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = LN.4(), type = "binomial")
  object0.W1.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = W1.4(), type = "binomial")
  object0.W2.4 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                      pmodels = list(~ treat - 1, ~ treat - 1,
                                     ~ 1, ~ treat -1),
                      data = data0, fct = W2.4(), type = "binomial")
  modelList0 <- list(object0.LL.4, object0.LN.4, object0.W1.4, object0.W2.4)
  
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  invisible(capture.output({
    resultKang <- suppressWarnings(bmdMA(modelList0, modelWeights = "Stacking", bmr = 0.77, def = "point", backgType = "modelBased", type = "Kang", stackingSeed = 1, stackingSplits = 3, display = FALSE))
    resultBoot <- suppressWarnings(bmdMA(modelList0, modelWeights = "Stacking", bmr = 0.77, def = "point", backgType = "modelBased", type = "bootstrap", R = 50, stackingSeed = 1, stackingSplits = 3, display = FALSE))
    }))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # resultKang
  expect_true(all(!is.na(resultKang$Results[, "BMD_MA"])))
  expect_equal(unname(resultKang$Results[, "BMD_MA"]), c(16.6498456454902, 30.3165644996408), tolerance = 1e-1)
  expect_equal(resultKang$Boot.samples.used, NA)
  expect_equal(unname(resultKang$interval[,"BMDL_MA"]), c(6.06814221070238,21.8988805204411), tolerance = 1)
  expect_equal(unname(resultKang$interval[,"BMDU_MA"]), c(27.231549080278,38.7342484788404), tolerance = 1)
  # resultBoot
  expect_true(all(!is.na(resultBoot$Results[, "BMD_MA"])))
  expect_equal(unname(resultBoot$Results[, "BMD_MA"]), c(16.6498456454902, 30.3165644996408), tolerance = 1e-1)
  expect_equal(resultBoot$Boot.samples.used, 50, tolerance = 1)
  expect_equal(unname(resultBoot$interval[,"BMDL_MA"]), c(2.12333205792999,25.08635184492), tolerance = 1)
  expect_equal(unname(resultBoot$interval[,"BMDU_MA"]), c(27.8840479507664,38.6902619548875), tolerance = 1)
})
