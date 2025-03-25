# Tests for MACurve function
# - Arguments and structure
#   - Missing arguments
# - aconiazide models


# Arguments and structure -------------------------------------------------

test_that("MACurve function handles missing required arguments", {
  aconiazide.LL.3 <- drm(weightChange ~ dose,data = drcData::aconiazide,fct = LL.3())
  aconiazide.W2.3 <- drm(weightChange ~ dose,data = drcData::aconiazide,fct = W2.3())
  
  expect_error(MACurve(modelList = list(aconiazide.LL.3, aconiazide.W2.3), modelWeights = "AIC"), 'argument "x" is missing, with no default')
  expect_error(MACurve(x = 1:10, modelWeights = "AIC"), 'argument "modelList" is missing, with no default')
  expect_error(MACurve(x = 1:10, modelList = list(aconiazide.LL.3, aconiazide.W2.3)), 'argument "modelWeights" is missing, with no default')
})


# aconiazide models -------------------------------------------------------

test_that("MACurve handles modelWeights argument", {
  aconiazide.LL.3 <- drm(weightChange ~ dose,data = drcData::aconiazide, fct = LL.3())
  aconiazide.LN.3 <- drm(weightChange ~ dose,data = drcData::aconiazide, fct = LN.3())
  aconiazide.W1.3 <- drm(weightChange ~ dose,data= drcData::aconiazide, fct = W1.3())
  aconiazide.W2.3 <- drm(weightChange ~ dose,data= drcData::aconiazide, fct = W2.3())
  modelList0 <- list(aconiazide.LL.3, aconiazide.LN.3,aconiazide.W1.3, aconiazide.W2.3)
  
  # model weights
  manWeights0 <- c(0.3, 0.2, 0.2, 0.3)
  AICWeights0 <- exp(-1/2 * (sapply(modelList0, AIC) - min(sapply(modelList0, AIC)))) / sum(exp(-1/2 * (sapply(modelList0, AIC) - min(sapply(modelList0, AIC)))) )
  BICWeights0 <- exp(-1/2 * sapply(modelList0, BIC)) / sum(exp(-1/2 * sapply(modelList0, BIC)) )
  set.seed(1)
  stackingWeights0 <- getStackingWeights(modelList0, nSplits = 3)
  
  x <- c(0, 50, 100, 150, 200, 250, 300, 350, 400)
  curveMat <- sapply(modelList0, function(mod) mod$curve[[1]](x))
  
  expect_equal(MACurve(x, modelList0, manWeights0), apply(curveMat, 1, function(z) sum(z * manWeights0)))
  expect_equal(MACurve(x, modelList0, "AIC"), apply(curveMat, 1, function(z) sum(z * AICWeights0)))
  expect_equal(MACurve(x, modelList0, "BIC"), apply(curveMat, 1, function(z) sum(z * BICWeights0)))
  expect_equal(MACurve(x, modelList0, "Stack", stackingSeed = 1, stackingSplits = 3), apply(curveMat, 1, function(z) sum(z * stackingWeights0)))
})
