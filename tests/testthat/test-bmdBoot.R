# Tests for bmdBoot function
# - Arguments and structure
#   - Missing arguments
#   - Correct backgType and def accepted
# - Simple model 
#   - correct bmd estimate 
# - Ryegrass model (continuous)
#   - correct bmd estimate (all definitions)
# - Ryegrass hormesis model (continuous)
#   - correct bmd estimate (all definitions)
# - TCDD model (binomial)
#   - correct bmd estimate (excess + additional)
# - Lemna model (count)
#   - correct bmd estimate (all definitions)
# - S.alba model (continuous with multiple curves)
#   - correct bmd estimate (point, extra, hybridExc)
# - Decreasing binomial model with multiple curves
#   - correct bmd estimate (point, extra, hybridExc)


# Arguments and structure -------------------------------------------------

test_that("bmdBoot function handles missing required arguments", {
  object0 <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2())
  
  expect_error(bmdBoot(), "object is missing")
  expect_error(bmdBoot(object0), "def is missing")
  expect_error(bmdBoot(lm(1:10 ~ 1)), 'object must be of class "drc"')
  expect_error(bmdBoot(object0, def = "invalid_def", backgType = "modelBased"), "Could not recognize def")
  expect_error(bmdBoot(object0, def = "excess", backgType = "invalid_type"), "Could not recognize backgType")
  expect_error(bmdBoot(object0, def = "excess"), "backgType is missing")
})


test_that("bmdBoot function accepts correct def", {
  object.cont <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2())
  object.binom <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2(), type = "binomial")
  object.poisson <- drm(y ~ x, data = data.frame(x = 1:5, y = c(12,11,3,0,0)), fct = LL.3(), type = "Poisson")
  
  # Binomial bmd def with continuous model
  expect_error(bmdBoot(object.cont, def = "excess", backgType = "modelBased", R = 1), '"excess" is not available for continuous data')
  expect_error(bmdBoot(object.cont, def = "additional", backgType = "modelBased", R = 1), '"additional" is not available for continuous data')
  
  # Binomial bmd def with Poisson model
  expect_error(bmdBoot(object.poisson, def = "excess", backgType = "modelBased", R = 1), '"excess" is not available for count data')
  expect_error(bmdBoot(object.poisson, def = "additional", backgType = "modelBased", R = 1), '"additional" is not available for count data')
  
  # Cont bmd def with binomial model
  expect_error(bmdBoot(object.binom, def = "relative", backgType = "modelBased", R = 1), '"relative" is not available for quantal data')
  expect_error(bmdBoot(object.binom, def = "extra", backgType = "modelBased", R = 1), '"extra" is not available for quantal data')
  expect_error(bmdBoot(object.binom, def = "added", backgType = "modelBased", R = 1), '"added" is not available for quantal data')
  expect_error(bmdBoot(object.binom, def = "hybridExc", backgType = "modelBased", R = 1), '"hybridExc" is not available for quantal data')
  expect_error(bmdBoot(object.binom, def = "hybridAdd", backgType = "modelBased", R = 1), '"hybridAdd" is not available for quantal data')
})

test_that("bmdBoot function returns expected structure", {
  object0 <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE, R = 1)
  
  expect_type(result, "list")
  expect_named(result, c("Results", "Boot.samples.used", "bootEst", "interval"))
  expect_s3_class(result, "bmd")
})


# Simple model results ----------------------------------------------------

test_that("bmdBoot function computes BMD (extra, bmr = 0.1) correctly for a simple model", {
  set.seed(1)
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 16.4578682695665)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(16.3490012773168,17.2565658001912))
})

test_that("bmdBoot function computes BMD (extra, bmr = 0.05) correctly for a simple model", {
  set.seed(1)
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmdBoot(object0, bmr = 0.05, def = "extra", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 15.4022544235763)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(15.2659263161485,16.4325341268289))
})

test_that("bmdBoot function computes BMD (extra, bmr = 0.1, bmdType = \"mean\") correctly for a simple model", {
  set.seed(1)
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE, R = 50, bmdType = "mean")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 16.780508604938)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(16.3490012773168,17.2565658001912))
})

test_that("bmdBoot function computes BMD (extra, bmr = 0.1, bmdType = \"median\") correctly for a simple model", {
  set.seed(1)
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE, R = 50, bmdType = "median")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 16.7232035648742)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(16.3490012773168,17.2565658001912))
})



# Ryegrass results --------------------------------------------------------

test_that("bmdBoot function computes BMD (point) correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 3.2, def = "point", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 3.64586140417992)
  expect_equal(unname(result$interval[1,]), c(3.26263598876729,3.96031450591884))
})

test_that("bmdBoot function computes BMD (extra) correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.46370565552042)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(1.34582379343345,1.65744656451026))
})

test_that("bmdBoot function computes BMD (relative) correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "relative", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.49902599632103)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(1.36986151499048,1.69766970621146))
})

test_that("bmdBoot function computes BMD (added) correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "added", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 0.728443033284576)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(0.601335381833136,0.938357650485433))
})

test_that("bmdBoot function computes BMD (hybridAdd with hybridSD background) correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.21255236145362)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(1.03230812002555,1.38259567213485))
})

test_that("bmdBoot function computes BMD (hybridExc with hybridSD background) correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.20672107998472)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(1.02698170147848,1.37683954154201))
})

test_that("bmdBoot function computes BMD (hybridExc with hybridPercentile background) correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.06888690340628)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(0.901388142148037,1.24291597301057))
})

test_that("bmdBoot function computes BMD (point, bootInterval = \"BCa\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 3.2, def = "point", backgType = "modelBased", bootInterval = "BCa", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 3.64586140417992)
  expect_equal(result$Results[1, "BMDL"], 3.40461643409508)
  expect_equal(result$interval[1,2], "Not available for BCa bootstrap")
})

test_that("bmdBoot function computes BMD (relative, bootInterval = \"BCa\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "relative", backgType = "modelBased", bootInterval = "BCa", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.49902599632103)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(result$Results[1,"BMDL"], 1.3364409073126)
  expect_equal(result$interval[1,2], "Not available for BCa bootstrap")
})

test_that("bmdBoot function computes BMD (hybridExc with hybridSD background, bootInterval = \"BCa\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, bootInterval = "BCa", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.20672107998472)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(result$Results[1,"BMDL"], 1.0046419139129)
  expect_equal(result$interval[1,2], "Not available for BCa bootstrap")
})

test_that("bmdBoot function computes BMD (hybridExc with hybridPercentile background, bootInterval = \"BCa\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, bootInterval = "BCa", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.06888690340628)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(result$Results[1,"BMDL"], 0.854123185248857)
  expect_equal(result$interval[1,2], "Not available for BCa bootstrap")
})

test_that("bmdBoot function computes BMD (point, bootType = \"parametric\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 3.2, def = "point", backgType = "modelBased", bootType = "parametric", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 3.64586140417992)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(3.12856944661618,4.34780971634149))
})

test_that("bmdBoot function computes BMD (relative, bootType = \"parametric\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "relative", backgType = "modelBased", bootType = "parametric", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.49902599632103)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(1.40708736919832,1.80364105456897))
})

test_that("bmdBoot function computes BMD (hybridExc with hybridSD background, bootType = \"parametric\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, bootType = "parametric", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.20672107998472)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(1.06044718186996,1.69830235608345))
})

test_that("bmdBoot function computes BMD (hybridExc with hybridPercentile background, bootType = \"parametric\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, bootType = "parametric", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.06888690340628)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(0.928960505221601,1.57848224390949))
})

test_that("bmdBoot function computes BMD (point, bootType = \"semiparametric\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 3.2, def = "point", backgType = "modelBased", bootType = "semiparametric", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 3.64586140417992)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(3.46643051130734,3.82733853806115))
})

test_that("bmdBoot function computes BMD (relative, bootType = \"semiparametric\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "relative", backgType = "modelBased", bootType = "semiparametric", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.49902599632103)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(1.33189055276109,1.70994880678072))
})

test_that("bmdBoot function computes BMD (hybridExc with hybridSD background, bootType = \"semiparametric\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, bootType = "semiparametric", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.20672107998472)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(0.942978153279434,1.36895504286135))
})

test_that("bmdBoot function computes BMD (hybridExc with hybridPercentile background, bootType = \"semiparametric\") correctly for ryegrass model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, bootType = "semiparametric", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.06888690340628)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(0.823230827522055,1.23451763714964))
})

test_that("bmdBoot function computes BMD (relative) correctly for ryegrass hormesis model", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "relative", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.55704870290614)
  expect_equal(result$Boot.samples.used, 42)
  expect_equal(unname(result$interval[1,]), c(1.44434658745943,1.75290180636397))
})

# test_that("bmdBoot function computes BMD (relative) with log-transformed response correctly for ryegrass model", {
#   set.seed(1)
#   object0 <- drm(log(rootl) ~ conc, data = drcData::ryegrass, fct = LL.4())
#   
#   result <- bmdBoot(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", display = FALSE, R = 50)
#   
#   # Expected results based on manual calculation (checked in v2.6.7)
#   expect_true(!is.na(result$Results[1, "BMD"]))
#   expect_equal(result$Results[1, "BMD"], 0.804218529940602)
#   expect_equal(result$Boot.samples.used, 0)
#   expect_equal(unname(result$interval[1,]), c(NA,NA))
# })
# 
# test_that("bmdBoot function computes BMD (relative) with square root-transformed response correctly for ryegrass model", {
#   set.seed(1)
#   object0 <- drm(sqrt(rootl) ~ conc, data = drcData::ryegrass, fct = LL.4())
#   
#   result <- bmdBoot(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", display = FALSE, R = 50)
#   
#   # Expected results based on manual calculation (checked in v2.6.7)
#   expect_true(!is.na(result$Results[1, "BMD"]))
#   expect_equal(result$Results[1, "BMD"], 1.29590294092622)
#   expect_equal(result$Boot.samples.used, 0)
#   expect_equal(unname(result$interval[1,]), c(NA,NA))
# })


test_that("bmdBoot function output remains consistent", {
  set.seed(1)
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  result <- bmdBoot(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE, R = 50)
  
  snapshot_data <- list(
    Results = as.list(result$Results),
    Boot.samples.used = as.list(result$Boot.samples.used),
    bootEst = as.list(result$bootEst),
    interval = as.list(result$interval)
  )
  
  # Store a snapshot of the entire result object
  expect_snapshot_value(snapshot_data, style = "deparse")
})



# TCDD results ------------------------------------------------------------

test_that("bmdBoot function computes BMD (point) correctly for TCDD model", {
  set.seed(1)
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- suppressWarnings(bmdBoot(object0, bmr = 0.22, def = "point", backgType = "modelBased", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 7.77184985530323)
  expect_equal(result$Boot.samples.used, 49)
  expect_equal(unname(result$interval[1,]), c(6.5722041232272,23.4873980148745))
})

test_that("bmdBoot function computes BMD (excess) correctly for TCDD model", {
  set.seed(1)
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- suppressWarnings(bmdBoot(object0, bmr = 0.05, def = "excess", backgType = "modelBased", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 5.56116921034511)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(1.28179749253742,4.68414395623765))
})

test_that("bmdBoot function computes BMD (additional) correctly for TCDD model", {
  set.seed(1)
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- suppressWarnings(bmdBoot(object0, bmr = 0.1, def = "additional", backgType = "modelBased", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 6.36475841679501)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(3.18111681473504,6.6921148746103))
})

test_that("bmdBoot function computes BMD (point, bootInterval = \"BCa\") correctly for TCDD model", {
  set.seed(1)
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- suppressWarnings(bmdBoot(object0, bmr = 0.22, def = "point", backgType = "modelBased", bootInterval = "BCa", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 7.77184985530323)
  expect_equal(result$Boot.samples.used, 49)
  expect_equal(result$Results[1, "BMDL"], 5.18330205847931)
  expect_equal(result$interval[1,2], "Not available for BCa bootstrap")
})

test_that("bmdBoot function computes BMD (excess, bootInterval = \"BCa\") correctly for TCDD model", {
  set.seed(1)
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- suppressWarnings(bmdBoot(object0, bmr = 0.05, def = "excess", backgType = "modelBased", bootInterval = "BCa", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 5.56116921034511)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(result$Results[1, "BMDL"], 6.14136910217948)
  expect_equal(result$interval[1,2], "Not available for BCa bootstrap")
})

test_that("bmdBoot function computes BMD (additional, bootInterval = \"BCa\") correctly for TCDD model", {
  set.seed(1)
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- suppressWarnings(bmdBoot(object0, bmr = 0.1, def = "additional", backgType = "modelBased", bootInterval = "BCa", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 6.36475841679501)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(result$Results[1, "BMDL"], 6.15359416425459)
  expect_equal(result$interval[1,2], "Not available for BCa bootstrap")
})

test_that("bmdBoot function computes BMD (point, bootType = \"parametric\") correctly for TCDD model", {
  set.seed(1)
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- suppressWarnings(bmdBoot(object0, bmr = 0.22, def = "point", backgType = "modelBased", bootType = "parametric", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 7.77184985530323)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(6.0472767939815,15.9893339617481))
})

test_that("bmdBoot function computes BMD (excess, bootType = \"parametric\") correctly for TCDD model", {
  set.seed(1)
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- suppressWarnings(bmdBoot(object0, bmr = 0.05, def = "excess", backgType = "modelBased", bootType = "parametric", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 5.56116921034511)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(1.28296661360405,4.74870734366396))
})

test_that("bmdBoot function computes BMD (additional, bootType = \"parametric\") correctly for TCDD model", {
  set.seed(1)
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- suppressWarnings(bmdBoot(object0, bmr = 0.1, def = "additional", backgType = "modelBased", bootType = "parametric", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 6.36475841679501)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(2.67704520143454,6.38843750211448))
})



# lemmna results ----------------------------------------------------------

test_that("bmdBoot function computes BMD (point) correctly for lemna model", {
  set.seed(1)
  object0 <- drm(frond.num ~ conc, data = drcData::lemna, fct = LL.3(), type = "Poisson")
  
  result <- bmdBoot(object0, bmr = 52, def = "point", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 4.35865965537475)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(3.6817392726952,4.96008552694743))
})

test_that("bmdBoot function computes BMD (extra) correctly for lemna model", {
  set.seed(1)
  object0 <- drm(frond.num ~ conc, data = drcData::lemna, fct = LL.3(), type = "Poisson")
  
  result <- bmdBoot(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 0.644966972651776)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(0.403304817944331,0.819989014246345))
})

test_that("bmdBoot function computes BMD (relative) correctly for lemna model", {
  set.seed(1)
  object0 <- drm(frond.num ~ conc, data = drcData::lemna, fct = LL.3(), type = "Poisson")
  
  result <- bmdBoot(object0, bmr = 0.1, def = "relative", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 0.644966972651776)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[1,]), c(0.403304817944331,0.819989014246345))
})


test_that("bmdBoot function computes BMD (point, bootInterval = \"BCa\") correctly for lemna model", {
  set.seed(1)
  object0 <- drm(frond.num ~ conc, data = drcData::lemna, fct = LL.3(), type = "Poisson")
  
  result <- bmdBoot(object0, bmr = 52, def = "point", backgType = "modelBased", bootInterval = "BCa", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 4.35865965537475)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(result$Results[1, "BMDL"], 3.89226934419482)
  expect_equal(result$interval[1,2], "Not available for BCa bootstrap")
})

test_that("bmdBoot function computes BMD (relative, bootInterval = \"BCa\") correctly for lemna model", {
  set.seed(1)
  object0 <- drm(frond.num ~ conc, data = drcData::lemna, fct = LL.3(), type = "Poisson")
  
  result <- bmdBoot(object0, bmr = 0.1, def = "relative", backgType = "modelBased", bootInterval = "BCa", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 0.644966972651776)
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(result$Results[1,"BMDL"], 0.46010142747602)
  expect_equal(result$interval[1,2], "Not available for BCa bootstrap")
})

test_that("bmdBoot function computes BMD (point, bootType = \"parametric\") correctly for lemna model", {
  set.seed(1)
  object0 <- drm(frond.num ~ conc, data = drcData::lemna, fct = LL.3(), type = "Poisson")
  
  expect_error(bmdBoot(object0, bmr = 52, def = "point", backgType = "modelBased", bootType = "parametric", display = FALSE, R = 50),
               "\"Poisson\" only works with nonparametric bootstrap")
})

test_that("bmdBoot function computes BMD (point, bootType = \"semiparametric\") correctly for lemna model", {
  set.seed(1)
  object0 <- drm(frond.num ~ conc, data = drcData::lemna, fct = LL.3(), type = "Poisson")
  
  expect_error(bmdBoot(object0, bmr = 52, def = "point", backgType = "modelBased", bootType = "semiparametric", display = FALSE, R = 50),
               "\"Poisson\" only works with nonparametric bootstrap")
})





# S.alba results ----------------------------------------------------------

test_that("bmdBoot function computes BMD (point) correctly for S.alba model", {
  set.seed(1)
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 3.2, def = "point", backgType = "modelBased", display = FALSE, R = 50)
  resultBCa <- suppressWarnings(bmdBoot(object0, bmr = 3.2, def = "point", backgType = "modelBased", bootInterval = "BCa", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(unname(result$Results[, "BMD"]), c(39.4912945056265, 22.1766859356908))
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[,"Lower"]), c(30.0422707998487,18.5473504361089))
  expect_equal(unname(result$interval[,"Upper"]), c(44.5736649922485,25.8729296047688))
  expect_true(all(!is.na(resultBCa$Results[, "BMD"])))
  expect_equal(unname(resultBCa$Results[, "BMD"]), c(39.4912945056265, 22.1766859356908))
  expect_equal(resultBCa$Boot.samples.used, 50)
  expect_equal(unname(resultBCa$Results[,"BMDL"]), c(34.0193865555056,19.45281610897))
  expect_equal(unname(resultBCa$interval[,"Upper"]), c("Not available for BCa bootstrap","Not available for BCa bootstrap"))
})

test_that("bmdBoot function computes BMD (point, bmdType = \"mean\") correctly for S.alba model", {
  set.seed(1)
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 3.2, def = "point", backgType = "modelBased", display = FALSE, R = 50, bmdType = "mean")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(unname(result$Results[, "BMD"]), c(37.5232019182286, 22.061872769014))
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[,"Lower"]), c(30.0422707998487,18.5473504361089))
  expect_equal(unname(result$interval[,"Upper"]), c(44.5736649922485,25.8729296047688))
})

test_that("bmdBoot function computes BMD (point, bmdType = \"median\") correctly for S.alba model", {
  set.seed(1)
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 3.2, def = "point", backgType = "modelBased", display = FALSE, R = 50, bmdType = "median")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(unname(result$Results[, "BMD"]), c(37.3570289171494, 21.8393124704904))
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[,"Lower"]), c(30.0422707998487,18.5473504361089))
  expect_equal(unname(result$interval[,"Upper"]), c(44.5736649922485,25.8729296047688))
})

test_that("bmdBoot function computes BMD (relative) correctly for S.alba model", {
  set.seed(1)
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.08, def = "relative", backgType = "modelBased", display = FALSE, R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(unname(result$Results[, "BMD"]), c(28.0790872125237, 18.9735396170819))
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[,"Lower"]), c(15.8093978884084,15.1505752899725))
  expect_equal(unname(result$interval[,"Upper"]), c(35.8713641261381,23.689950459446))
})

test_that("bmdBoot function computes BMD (hybridExc with hybridSD background) correctly for S.alba model", {
  set.seed(1)
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  
  result <- bmdBoot(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, display = FALSE, R = 50)
  resultBCa <- suppressWarnings(bmdBoot(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, bootInterval = "BCa", display = FALSE, R = 50))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(unname(result$Results[, "BMD"]), c(28.0253530227688, 19.0291591246355))
  expect_equal(result$Boot.samples.used, 50)
  expect_equal(unname(result$interval[,"Lower"]), c(13.5226993261253,14.6581794833924))
  expect_equal(unname(result$interval[,"Upper"]), c(36.5969418859266,22.9059472014908))
  expect_true(all(!is.na(resultBCa$Results[, "BMD"])))
  # resultBCa
  expect_equal(unname(resultBCa$Results[, "BMD"]), c(28.0253530227688, 19.0291591246355))
  expect_equal(resultBCa$Boot.samples.used, 50)
  expect_equal(unname(resultBCa$Results[,"BMDL"]), c(17.8777299232804,14.558986489691))
  expect_equal(unname(resultBCa$interval[,"Upper"]), c("Not available for BCa bootstrap","Not available for BCa bootstrap"))
})

test_that("bmdBoot function output remains consistent with model with multiple curves", {
  set.seed(1)
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  result <- bmdBoot(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE, R = 50)
  
  snapshot_data <- list(
    Results = as.list(result$Results),
    Boot.samples.used = as.list(result$Boot.samples.used),
    bootEst = as.list(result$bootEst),
    interval = as.list(result$interval)
  )
  
  # Store a snapshot of the entire result object
  expect_snapshot_value(snapshot_data, style = "deparse")
})


# Decreasing binomial model with multiple curves --------------------------
test_that("bmdBoot function computes BMD (point) correctly for TCDD model", {
  data0 <- data.frame(
    conc = c(0, rep(c(20, 40, 80, 160, 320), 2)),
    total = c(50, rep(20, 5*2)),
    alive = c(47, 14, 11, 6, 9, 6, 19, 11, 8, 5, 3),
    treat = c("Control", rep(c("Treat1", "Treat2"), each = 5))
  )
  
  object0 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                 pmodels = list(~ treat - 1, ~ treat - 1,
                                ~ 1, ~ treat -1),
                 data = data0, fct = W2.4(), type = "binomial")
  
  set.seed(1)
  invisible(capture.output({
    result <- suppressWarnings(bmdBoot(object0, bmr = 0.77, def = "point", backgType = "modelBased", display = FALSE, R = 50))
    resultBCa <- suppressWarnings(bmdBoot(object0, bmr = 0.77, def = "point", backgType = "modelBased", bootInterval = "BCa", display = FALSE, R = 50))
  }))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(unname(result$Results[, "BMD"]), c(16.7094280948318, 30.2464214935383))
  expect_equal(result$Boot.samples.used, 44)
  expect_equal(unname(result$interval[,"Lower"]), c(2.80478996663512,25.4932116970026))
  expect_equal(unname(result$interval[,"Upper"]), c(20.9736682163773,35.0170372426051))
  # resultBCa
  expect_true(all(!is.na(resultBCa$Results[, "BMD"])))
  expect_equal(unname(resultBCa$Results[, "BMD"]), c(16.7094280948318, 30.2464214935383))
  expect_equal(resultBCa$Boot.samples.used, 44)
  expect_equal(unname(resultBCa$Results[,"BMDL"]), c(12.2322036611872,25.3303655120452))
  expect_equal(unname(resultBCa$interval[,"Upper"]), c("Not available for BCa bootstrap","Not available for BCa bootstrap"))
})

test_that("bmdBoot function computes BMD (excess) correctly for TCDD model", {
  data0 <- data.frame(
    conc = c(0, rep(c(20, 40, 80, 160, 320), 2)),
    total = c(50, rep(20, 5*2)),
    alive = c(47, 14, 11, 6, 9, 6, 19, 11, 8, 5, 3),
    treat = c("Control", rep(c("Treat1", "Treat2"), each = 5))
  )
  
  object0 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                 pmodels = list(~ treat - 1, ~ treat - 1,
                                ~ 1, ~ treat -1),
                 data = data0, fct = W2.4(), type = "binomial")
  
  set.seed(1)
  invisible(capture.output({
    result <- suppressWarnings(bmdBoot(object0, bmr = 0.1, def = "excess", backgType = "modelBased", display = FALSE, R = 50))
    resultBCa <- suppressWarnings(bmdBoot(object0, bmr = 0.1, def = "excess", backgType = "modelBased", bootInterval = "BCa", display = FALSE, R = 50))
  }))
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(unname(result$Results[, "BMD"]), c(12.7945107782873, 25.107888955923))
  expect_equal(result$Boot.samples.used, 45)
  expect_equal(unname(result$interval[,"Lower"]), c(0.0803383339260978,15.8070091895054))
  expect_equal(unname(result$interval[,"Upper"]), c(16.4173994133493,29.5166318382001))
  # resultBCa
  expect_true(all(!is.na(resultBCa$Results[, "BMD"])))
  expect_equal(unname(resultBCa$Results[, "BMD"]), c(12.7945107782873, 25.107888955923))
  expect_equal(resultBCa$Boot.samples.used, 46)
  expect_equal(unname(resultBCa$Results[,"BMDL"]), c(14.1406750203466,17.9700937734153))
  expect_equal(unname(resultBCa$interval[,"Upper"]), c("Not available for BCa bootstrap","Not available for BCa bootstrap"))
})

test_that("bmdBoot function computes BMD (additional) correctly for TCDD model", {
  data0 <- data.frame(
    conc = c(0, rep(c(20, 40, 80, 160, 320), 2)),
    total = c(50, rep(20, 5*2)),
    alive = c(47, 14, 11, 6, 9, 6, 19, 11, 8, 5, 3),
    treat = c("Control", rep(c("Treat1", "Treat2"), each = 5))
  )
  
  object0 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                 pmodels = list(~ treat - 1, ~ treat - 1,
                                ~ 1, ~ treat -1),
                 data = data0, fct = W2.4(), type = "binomial")
  
  set.seed(1)
  invisible(capture.output(
    result <- suppressWarnings(bmdBoot(object0, bmr = 0.1, def = "additional", backgType = "modelBased", display = FALSE, R = 50))
  ))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(unname(result$Results[, "BMD"]), c(13.0508853789932, 25.4655976138905))
  expect_equal(result$Boot.samples.used, 45)
  expect_equal(unname(result$interval[,"Lower"]), c(0.0542035268237113,15.6266599428455))
  expect_equal(unname(result$interval[,"Upper"]), c(16.8396337661858,30.3505265877131))
})
