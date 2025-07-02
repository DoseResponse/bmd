# Tests for bmdHetVar function
# - Arguments and structure
#   - Missing arguments
#   - Correct backgType and def accepted
# - Ryegrass model
#   - correct bmd estimate (hybridExc + hybridAdd with hybridPercentile and hybridSD backgType)
# - GiantKelp model
#   - correct bmd estimate (hybridExc + hybridAdd with hybridPercentile and hybridSD backgType)




# Arguments and structure -------------------------------------------------

test_that("bmdHetVarMA handles missing required arguments", {
  lm_object <- lm(y ~ x, 
                  data = data.frame(x = 0:4, 
                                    y = 1:5 + c(-0.4, 0, 0.2, -0.1, 0.13)))
  ryegrass.List <- list(drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4()),
                        drm(rootl ~ conc, data = drcData::ryegrass, fct = W2.4()))
  var.formula0 <- ~ fitted + I(fitted^2)
  ryegrass.hetVarList <- list(drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = LL.4()),
                              drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W2.4()))
  
  expect_error(bmdHetVarMA(lm_object), 'modelList must be a list of dose-response models with a heterogeneous variance structure of class "drcHetVar"')
  expect_error(bmdHetVarMA(ryegrass.List), 'modelList must be a list of dose-response models with a heterogeneous variance structure of class "drcHetVar"')
  expect_error(bmdHetVarMA(ryegrass.hetVarList), 'modelWeights must either be "AIC", "BIC" or a numeric vector of same length as modelList')
  expect_error(bmdHetVarMA(ryegrass.hetVarList, modelWeights = "noWeights"), 'modelWeights must either be "AIC", "BIC" or a numeric vector of same length as modelList')
  expect_error(bmdHetVarMA(ryegrass.hetVarList, modelWeights = c(1,2,3)), 'modelWeights must either be "AIC", "BIC" or a numeric vector of same length as modelList')
  expect_error(bmdHetVarMA(ryegrass.hetVarList, modelWeights = "AIC", bmr = 0.1, def = "hybridExc"), 'backgType is missing. Options are "absolute", "hybridSD" or "hybridPercentile"')
  expect_error(bmdHetVarMA(ryegrass.hetVarList, modelWeights = "AIC", bmr = 0.1, def = "hybridExc", backgType = "wrongType"), 'Could not recognize backgType. Options are "absolute", "hybridSD" or "hybridPercentile"')
  expect_error(bmdHetVarMA(bmr = 0.1, backgType = "hybridPercentile", def = "hybridExc"), 'argument "modelList" is missing, with no default')
  expect_error(bmdHetVarMA(ryegrass.hetVarList, modelWeights = "AIC", backgType = "hybridPercentile", def = "hybridExc"), 'argument "bmr" needs to be specified as a number between 0 and 1')
  expect_error(bmdHetVarMA(ryegrass.hetVarList, modelWeights = "AIC", bmr = 1.5, backgType = "hybridPercentile", def = "hybridExc"), 'argument "bmr" needs to be specified as a number between 0 and 1')
  expect_error(bmdHetVarMA(ryegrass.hetVarList, modelWeights = "AIC", bmr = 0.1, backgType = "hybridPercentile"), 'def is missing. Options are "hybridExc" or "hybridAdd"')
  expect_error(bmdHetVarMA(ryegrass.hetVarList, modelWeights = "AIC", bmr = 0.1, backgType = "hybridPercentile", def = "wrongDef"), 'Could not recognize def. Options are "hybridExc" or "hybridAdd"')
})





# Ryegrass model ----------------------------------------------------------

test_that("bmdHetVarMA on Ryegrass models, def = hybridExc, backgtype = hybridPercentile", {
  var.formula0 <- ~ fitted + I(fitted^2)
  modelList0 <- list(drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = LL.4()),
                     drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = LN.4()),
                     drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W1.4()),
                     drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W2.4()))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultAIC <- bmdHetVarMA(modelList0, modelWeights = "AIC", bmr = 0.1, backgType = "hybridPercentile", 
                        backg = 0.1, def = "hybridExc", interval = "boot", R = 50, level = 0.95, 
                        progressInfo = FALSE, display = FALSE)
  resultBIC <- bmdHetVarMA(modelList0, modelWeights = "BIC", bmr = 0.1, backgType = "hybridPercentile", 
                           backg = 0.1, def = "hybridExc", interval = "boot", R = 50, level = 0.95, 
                           progressInfo = FALSE, display = FALSE)
  resultManWeights <- bmdHetVarMA(modelList0, modelWeights = c(0.2, 0.1, 0.1, 0.6), bmr = 0.1, backgType = "hybridPercentile", 
                                  backg = 0.1, def = "hybridExc", interval = "boot", R = 50, level = 0.95, 
                                  progressInfo = FALSE, display = FALSE)
  
  # resultAIC
  expect_true(!is.na(resultAIC$Results[1, "BMD_MA"]))
  expect_equal(resultAIC$Results[1, "BMD_MA"], c(BMD_MA=1.09399583860305))
  expect_equal(resultAIC$Boot.samples.used, 48)
  expect_equal(unname(resultAIC$interval[1,]), c(0.91243677103485,1.3064583554297), tolerance = 1e-1)
  expect_equal(resultAIC$modelWeights, c(0.239817917878516,0.112398486139835,0.00397069684089491,0.643812899140755), tolerance = 1e-4)
  
  # resultBIC
  expect_true(!is.na(resultBIC$Results[1, "BMD_MA"]))
  expect_equal(resultBIC$Results[1, "BMD_MA"], c(BMD_MA=1.09399583860305))
  expect_equal(resultBIC$Boot.samples.used, 49)
  expect_equal(unname(resultBIC$interval[1,]), c(0.740131127613234,1.29617185295069), tolerance = 1e-1)
  expect_equal(resultBIC$modelWeights, c(0.239817917878515,0.112398486139835,0.00397069684089491,0.643812899140755), tolerance = 1e-4)
  
  # resultManWeights
  expect_true(!is.na(resultManWeights$Results[1, "BMD_MA"]))
  expect_equal(resultManWeights$Results[1, "BMD_MA"], c(BMD_MA=1.07129213615452))
  expect_equal(resultManWeights$Boot.samples.used, 50)
  expect_equal(unname(resultManWeights$interval[1,]), c(0.95908147844657,1.27941263725323), tolerance = 1e-1)
  expect_equal(resultManWeights$modelWeights, c(0.2,0.1,0.1,0.6))
})

test_that("bmdHetVarMA on Ryegrass models, def = hybridExc, backgtype = hybridSD", {
  var.formula0 <- ~ fitted + I(fitted^2)
  modelList0 <- list(drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = LL.4()),
                  drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W2.4()))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVarMA(modelList0, modelWeights = c(0.4, 0.6), bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridExc", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD_MA"]))
  expect_equal(result$Results[1, "BMD_MA"], c(BMD_MA=1.2904405178104))
  expect_equal(result$Boot.samples.used, 49)
  expect_equal(unname(result$interval[1,]), c(1.1544076635982,1.5029430687596), tolerance = 1e-1)
  expect_equal(result$modelWeights, c(0.4,0.6))
})

test_that("bmdHetVarMA on Ryegrass models, def = hybridAdd, backgtype = hybridPercentile", {
  modelList0 <- list(drmHetVar(rootl ~ conc, ~ fitted + I(fitted^2), data = drcData::ryegrass, fct = LL.4()),
                     drmHetVar(rootl ~ conc, ~ fitted, data = drcData::ryegrass, fct = LL.4()))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVarMA(modelList0, modelWeights = "AIC", bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD_MA"]))
  expect_equal(result$Results[1, "BMD_MA"], c(BMD_MA=0.856751513601129))
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(0.662428472523196,1.10812296172592), tolerance = 1e-1)
  expect_equal(result$modelWeights, c(0.979658442072015,0.0203415579279849))
})

test_that("bmdHetVarMA on Ryegrass models, def = hybridAdd, backgtype = hybridSD", {
  var.formula0 <- ~ fitted + I(fitted^2)
  modelList0 <- list(drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = LL.4()),
                     drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W2.4()))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVarMA(modelList0, modelWeights = "AIC", bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD_MA"]))
  expect_equal(result$Results[1, "BMD_MA"], c(BMD_MA=1.33134958361854))
  expect_equal(result$Boot.samples.used, 49)
  expect_equal(unname(result$interval[1,]), c(1.14730351933847,1.51123369350389), tolerance = 1e-1)
  expect_equal(result$modelWeights, c(0.271400581848749,0.728599418151251))
})



# GiantKelp model ---------------------------------------------------------


test_that("bmdHetVarMA on GiantKelp models, def = hybridExc, backgtype = hybridPercentile", {
  modelList0 <- list(drmHetVar(tubeLength ~ dose, ~ log(dose+1) + I(log(dose+1)^2), data = drcData::GiantKelp, fct = LL.4()),
                     drmHetVar(tubeLength ~ dose, ~ fitted + I(fitted^2), data = drcData::GiantKelp, fct = LL.4()))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  resultAIC <- bmdHetVarMA(modelList0, modelWeights = "AIC", bmr = 0.1, backgType = "hybridPercentile", 
                           backg = 0.1, def = "hybridExc", interval = "boot", R = 50, level = 0.95, 
                           progressInfo = FALSE, display = FALSE)
  resultBIC <- bmdHetVarMA(modelList0, modelWeights = "BIC", bmr = 0.1, backgType = "hybridPercentile", 
                           backg = 0.1, def = "hybridExc", interval = "boot", R = 50, level = 0.95, 
                           progressInfo = FALSE, display = FALSE)
  resultManWeights <- bmdHetVarMA(modelList0, modelWeights = c(0.3, 0.7), bmr = 0.1, backgType = "hybridPercentile", 
                                  backg = 0.1, def = "hybridExc", interval = "boot", R = 50, level = 0.95, 
                                  progressInfo = FALSE, display = FALSE)
  # resultAIC
  expect_true(!is.na(resultAIC$Results[1, "BMD_MA"]))
  expect_equal(resultAIC$Results[1, "BMD_MA"], c(BMD_MA=5.22750504984351))
  expect_equal(resultAIC$Boot.samples.used, 48)
  expect_equal(unname(resultAIC$interval[1,]), c(0.970896171135017,17.9687826158524), tolerance = 1e-1)
  expect_equal(resultAIC$modelWeights, c(0.75212289732896,0.24787710267104))
  
  # resultBIC
  expect_true(!is.na(resultBIC$Results[1, "BMD_MA"]))
  expect_equal(resultBIC$Results[1, "BMD_MA"], c(BMD_MA=5.22750504984351))
  expect_equal(resultBIC$Boot.samples.used, 50)
  expect_equal(unname(resultBIC$interval[1,]), c(1.93813846559323,15.0689738215226), tolerance = 1e-1)
  expect_equal(resultBIC$modelWeights, c(0.75212289732896,0.24787710267104))
  
  # resultManWeights
  expect_true(!is.na(resultManWeights$Results[1, "BMD_MA"]))
  expect_equal(resultManWeights$Results[1, "BMD_MA"], c(BMD_MA=4.49558070194896))
  expect_equal(resultManWeights$Boot.samples.used, 50)
  expect_equal(unname(resultManWeights$interval[1,]), c(0.509778477285164,16.6029584852476), tolerance = 1e-1)
  expect_equal(resultManWeights$modelWeights, c(0.3,0.7))
})

test_that("bmdHetVarMA on GiantKelp models, def = hybridExc, backgtype = hybridSD", {
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  modelList0 <- list(drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LL.4()),
                     drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LN.4()))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVarMA(modelList0, modelWeights = c(0.6, 0.4), bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridExc", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD_MA"]))
  expect_equal(result$Results[1, "BMD_MA"], c(BMD_MA=10.4673358308111))
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(3.46365031231543,20.4834802123576), tolerance = 1e-1)
  expect_equal(result$modelWeights, c(0.6,0.4))
})

test_that("bmdHetVarMA on GiantKelp models, def = hybridAdd, backgtype = hybridPercentile", {
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  modelList0 <- list(drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LL.4()),
                     drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LN.4()))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVarMA(modelList0, modelWeights = "AIC", bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD_MA"]))
  expect_equal(result$Results[1, "BMD_MA"], c(BMD_MA=6.3975810199722))
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(1.48403844387169,17.2022294707629), tolerance = 1e-1)
  expect_equal(result$modelWeights, c(0.5788385115668,0.421161488433199))
})

test_that("bmdHetVarMA on GiantKelp models, def = hybridAdd, backgtype = hybridSD", {
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  modelList0 <- list(drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LL.4()),
                     drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LN.4()))
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVarMA(modelList0, modelWeights = c(0.6, 0.4), bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD_MA"]))
  expect_equal(result$Results[1, "BMD_MA"], c(BMD_MA=10.5823318236905))
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(3.55888343309975,20.6224431910256), tolerance = 1e-1)
  expect_equal(result$modelWeights, c(0.6,0.4))
})



