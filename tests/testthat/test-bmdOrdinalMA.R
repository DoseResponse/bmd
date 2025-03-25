# Tests for bmdOrdinalMA function
# - Arguments and structure
#   - Missing arguments
#   - Correct backgType and def accepted
# - guthionS model
#   - correct bmd estimate (all definitions)
#   - delta and bootstrap intervals




# Arguments and structure -------------------------------------------------

## TO BE ADDED




# guthionS model ----------------------------------------------------------

test_that("bmdOrdinalMA computes correct bmd estimates (def = point) for gutionS model", {
  guthionS <- subset(drcData::guthion, trt == "S")

  guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = LL.2())
  guthionS.W1 <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = W1.2())
  guthionS.W2 <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = W2.2())
  modelList0 <- list(guthionS.LL, guthionS.W1, guthionS.W2)
  manWeights0 <- c(0.3, 0.6, 0.1)
  resultKangAIC <- bmdOrdinalMA(modelList0, modelWeights = "AIC", bmr = 0.2, backgType = "modelBased", def = "point", type = "Kang", display = FALSE)
  resultKangBIC <- bmdOrdinalMA(modelList0, modelWeights = "BIC", bmr = 0.2, backgType = "modelBased", def = "point", type = "Kang", display = FALSE)
  resultKangManWeights <- bmdOrdinalMA(modelList0, modelWeights = manWeights0, bmr = 0.2, backgType = "modelBased", def = "point", type = "Kang", display = FALSE)
  set.seed(1)
  resultBootAIC <- bmdOrdinalMA(modelList0, modelWeights = "AIC", bmr = 0.2, backgType = "modelBased", def = "point", type = "bootstrap", R = 50, progressInfo = FALSE, display = FALSE)
  resultBootBIC <- bmdOrdinalMA(modelList0, modelWeights = "BIC", bmr = 0.2, backgType = "modelBased", def = "point", type = "bootstrap", R = 50, progressInfo = FALSE, display = FALSE)
  resultBootManWeights <- bmdOrdinalMA(modelList0, modelWeights = manWeights0, bmr = 0.2, backgType = "modelBased", def = "point", type = "bootstrap", R = 50, progressInfo = FALSE, display = FALSE)

  # checks
  # resultKangAIC
  expect_true(all(!is.na(resultKangAIC$Results[,1])))
  expect_equal(resultKangAIC$Results[,1], c("moribund+dead" = 23.2544468633164, "dead" = 24.3204043959293))
  expect_equal(resultKangAIC$Results[,2], c("moribund+dead" = 20.9889153780884, "dead" = 21.6671551612867))
  expect_equal(resultKangAIC$interval[,2], c("moribund+dead" = 25.5199783485443, "dead" = 26.9736536305719))
  expect_true(is.na(resultKangAIC$Boot.samples.used))
  # resultKangBIC
  expect_true(all(!is.na(resultKangBIC$Results[,1])))
  expect_equal(resultKangBIC$Results[,1], c("moribund+dead" = 23.2544468633164, "dead" = 24.3204043959293))
  expect_equal(resultKangBIC$Results[,2], c("moribund+dead" = 20.9889153780884, "dead" = 21.6671551612867))
  expect_equal(resultKangBIC$interval[,2], c("moribund+dead" = 25.5199783485443, "dead" = 26.9736536305719))
  expect_true(is.na(resultKangBIC$Boot.samples.used))
  # resultKangManWeights
  expect_true(all(!is.na(resultKangManWeights$Results[,1])))
  expect_equal(resultKangManWeights$Results[,1], c("moribund+dead" = 23.7749220575841, "dead" = 24.7506554039273))
  expect_equal(resultKangManWeights$Results[,2], c("moribund+dead" = 21.2467089545517, "dead" = 21.830241980314))
  expect_equal(resultKangManWeights$interval[,2], c("moribund+dead" = 26.3031351606164, "dead" = 27.6710688275405))
  expect_true(is.na(resultKangManWeights$Boot.samples.used))
  # resultBootAIC
  expect_equal(resultBootAIC$Results[,1], c("moribund+dead" = 23.2544468633164, "dead" = 24.3204043959293))
  expect_equal(resultBootAIC$Results[,2], c("moribund+dead" = 22.25059194107, "dead" = 22.8512766871584))
  expect_equal(resultBootAIC$interval[,2], c("moribund+dead" = 28.6789792473312, "dead" = 29.1978148632589))
  expect_equal(resultBootAIC$Boot.samples.used, 50)
  # resultBootBIC
  expect_equal(resultBootBIC$Results[,1], c("moribund+dead" = 23.2544468633164, "dead" = 24.3204043959293))
  expect_equal(resultBootBIC$Results[,2], c("moribund+dead" = 21.3030492748921, "dead" = 22.8711536934761))
  expect_equal(resultBootBIC$interval[,2], c("moribund+dead" = 28.2940621821775, "dead" = 29.467583602007))
  expect_equal(resultBootBIC$Boot.samples.used, 50)
  # resultBootManWeights
  expect_equal(resultBootManWeights$Results[,1], c("moribund+dead" = 23.7749220575841, "dead" = 24.7506554039273))
  expect_equal(resultBootManWeights$Results[,2], c("moribund+dead" = 21.9772468661395, "dead" = 21.8972828765798))
  expect_equal(resultBootManWeights$interval[,2], c("moribund+dead" = 26.1263977743362, "dead" = 27.4023387249217))
  expect_equal(resultBootManWeights$Boot.samples.used, 50)
})

test_that("bmdOrdinal computes correct bmd estimates (def = excess) for gutionS model", {
  guthionS <- subset(drcData::guthion, trt == "S")

  guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = LL.2())
  guthionS.W1 <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = W1.2())
  guthionS.W2 <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = W2.2())
  modelList0 <- list(guthionS.LL, guthionS.W1, guthionS.W2)
  
  result <- bmdOrdinalMA(modelList0, modelWeights = "AIC", bmr = 0.1, backgType = "modelBased", def = "excess", type = "Kang", display = FALSE)

  # checks
  expect_true(all(!is.na(result$Results[,1])))
  expect_equal(result$Results[,1], c("moribund+dead" = 20.0639900692681, "dead" = 20.4469882115098))
  expect_equal(result$Results[,2], c("moribund+dead" = 17.6456398082174, "dead" = 17.5754769224672))
  expect_equal(result$interval[,2], c("moribund+dead" = 22.4823403303188, "dead" = 23.3184995005524))
})

test_that("bmdOrdinal computes correct bmd estimates (def = additional) for gutionS model", {
  guthionS <- subset(drcData::guthion, trt == "S")
  
  guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = LL.2())
  guthionS.W1 <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = W1.2())
  guthionS.W2 <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = W2.2())
  modelList0 <- list(guthionS.LL, guthionS.W1, guthionS.W2)
  
  result <- bmdOrdinalMA(modelList0, modelWeights = "AIC", bmr = 0.1, backgType = "modelBased", def = "additional", type = "Kang", display = FALSE)

  # checks
  expect_true(all(!is.na(result$Results[,1])))
  expect_equal(result$Results[,1], c("moribund+dead" = 20.0639900692681, "dead" = 20.4469882115098))
  expect_equal(result$Results[,2], c("moribund+dead" = 17.6456398082174, "dead" = 17.5754769224672))
  expect_equal(result$interval[,2], c("moribund+dead" = 22.4823403303188, "dead" = 23.3184995005524))
})


