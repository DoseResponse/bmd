# Tests for bmd function
# - Arguments and structure
#   - Missing arguments
#   - Correct backgType and def accepted
# - Simple model 
#   - correct bmd estimate 
#   - delta and profile intervals
# - Ryegrass model (continuous)
#   - correct bmd estimate (all definitions)
#   - delta, inv and profile intervals
# - Ryegrass hormesis model (continuous)
#   - correct bmd estimate (all definitions)
#   - delta and profile intervals
# - TCDD model (binomial)
#   - correct bmd estimate (excess + additional)
#   - delta and inv intervals
# - Chlorac model (binomial)
#   - correct bmd estimate (excess + additional)
#   - delta, profile and inv intervals
# - Lemna model (count)
#   - correct bmd estimate (all definitions)
#   - delta, inv and profile intervals
# - S.alba model (continuous with multiple curves)
#   - correct bmd estimate (point, extra, hybridExc)
#   - delta and profile intervals
# - Increasing continuous model
#   - correct bmd estimate (all definitions)
#   - delta, inv and profile intervals
# - Decreasing binomial model with multiple curves
#   - correct bmd estimate (point, extra, hybridExc)
#   - delta
# - Meta analytic model (drmMMRE)


# Arguments and structure -------------------------------------------------

test_that("bmd function handles missing required arguments", {
  object0 <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2())
  
  expect_error(bmd(), "object is missing")
  expect_error(bmd(object0), "def is missing")
  expect_error(bmd(lm(1:10 ~ 1)), 'object must be of class "drc"')
  expect_error(bmd(object0, def = "invalid_def", backgType = "modelBased"), "Could not recognize def")
  expect_error(bmd(object0, def = "excess", backgType = "invalid_type"), "Could not recognize backgType")
  expect_error(bmd(object0, def = "excess"), "\"excess\" is not available for continuous data")
})


test_that("bmd function accepts correct def", {
  object.cont <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2())
  object.binom <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2(), type = "binomial")
  object.poisson <- drm(y ~ x, data = data.frame(x = 1:5, y = c(12,11,3,0,0)), fct = LL.3(), type = "Poisson")
  
  # Binomial bmd def with continuous model
  expect_error(bmd(object.cont, def = "excess", backgType = "modelBased"), '"excess" is not available for continuous data')
  expect_error(bmd(object.cont, def = "additional", backgType = "modelBased"), '"additional" is not available for continuous data')
  
  # Binomial bmd def with Poisson model
  expect_error(bmd(object.poisson, def = "excess", backgType = "modelBased"), '"excess" is not available for count data')
  expect_error(bmd(object.poisson, def = "additional", backgType = "modelBased"), '"additional" is not available for count data')
  
  # Cont bmd def with binomial model
  expect_error(bmd(object.binom, def = "relative", backgType = "modelBased"), '"relative" is not available for quantal data')
  expect_error(bmd(object.binom, def = "extra", backgType = "modelBased"), '"extra" is not available for quantal data')
  expect_error(bmd(object.binom, def = "added", backgType = "modelBased"), '"added" is not available for quantal data')
  expect_error(bmd(object.binom, def = "hybridExc", backgType = "modelBased"), '"hybridExc" is not available for quantal data')
  expect_error(bmd(object.binom, def = "hybridAdd", backgType = "modelBased"), '"hybridAdd" is not available for quantal data')
})

test_that("bmd function returns expected structure", {
  object0 <- drm(y ~ x, data = data.frame(x = 1:5, y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  
  expect_type(result, "list")
  expect_named(result, c("Results", "bmrScaled", "interval", "SE", "model"))
  expect_s3_class(result, "bmd")
})


# Simple model results ----------------------------------------------------

test_that("bmd function computes BMD (extra, bmr = 0.1) correctly for a simple model", {
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 16.4578682695665)
  expect_equal(result$bmrScaled[1], 0.9)
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
})

test_that("bmd function computes BMD (extra, bmr = 0.05) correctly for a simple model", {
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmd(object0, bmr = 0.05, def = "extra", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 15.4022544235763)
  expect_equal(result$bmrScaled[1], 0.95)
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
})

test_that("bmd function computes correct confidence interval for a simple model", {
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$interval[1, "Lower"]))
  expect_true(!is.na(result$interval[1, "Upper"]))
  expect_equal(result$Results[1, "BMDL"], result$interval[1, "Lower"])
  expect_equal(result$interval[1, "Lower"], 16.0631524054276)
  expect_equal(result$interval[1, "Upper"], 16.8525841337055)
})

test_that("bmd function computes correct profile confidence interval for a simple model", {
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$interval[1, "Lower"]))
  expect_true(!is.na(result$interval[1, "Upper"]))
  expect_equal(result$Results[1, "BMDL"], result$interval[1, "Lower"])
  expect_equal(result$interval[1, "Lower"], 16.3431132583152)
  expect_equal(result$interval[1, "Upper"], 17.2291662734977)
})

test_that("bmd function computes correct profile grid confidence interval for a simple model", {
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "profileGrid", profileGridSize = 50, profileProgressInfo = FALSE, display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$interval[1, "Lower"]))
  expect_true(!is.na(result$interval[1, "Upper"]))
  expect_equal(result$Results[1, "BMDL"], result$interval[1, "Lower"])
  expect_equal(result$interval[1, "Lower"], 16.3452825544063)
  expect_equal(result$interval[1, "Upper"], 16.968912613503)
})

test_that("bmd function computes correct inverse regression confidence interval for a simple model", {
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "inv", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$interval[1, "Lower"]))
  expect_true(!is.na(result$interval[1, "Upper"]))
  expect_equal(result$Results[1, "BMDL"], result$interval[1, "Lower"])
  expect_equal(result$interval[1, "Lower"], 16.0671734467641)
  expect_equal(result$interval[1, "Upper"], 16.8460197017419)
})

test_that("bmd function computes correct sandwich confidence interval for a simple model", {
  object0 <- drm(y ~ x, data = data.frame(x = c(0,10,20,40,80), y = c(1,1,0.5,0,0)), fct = LL.2())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "sandwich", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$interval[1, "Lower"]))
  expect_true(!is.na(result$interval[1, "Upper"]))
  expect_equal(result$Results[1, "BMDL"], result$interval[1, "Lower"])
  expect_equal(result$interval[1, "Lower"], 16.1160251608648)
  expect_equal(result$interval[1, "Upper"], 16.7997113782683)
})


# Ryegrass results --------------------------------------------------------

test_that("bmd function computes BMD (point) correctly for ryegrass model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmd(object0, bmr = 3.2, def = "point", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 3.2, def = "point", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultProfile <- bmd(object0, bmr = 3.2, def = "point", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 3.64586140417992)
  expect_equal(result$bmrScaled[1,1], 3.2)
  expect_equal(unname(result$interval[1,]), c(3.29979318251471,3.99192962584513))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  # sandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 3.64586140417992)
  expect_equal(resultSandwich$bmrScaled[1,1], 3.2)
  expect_equal(unname(resultSandwich$interval[1,]), c(3.05041199694892,4.24131081141093))
  # profile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 3.64586140417992)
  expect_equal(resultProfile$bmrScaled[1,1], 3.2)
  expect_equal(unname(resultProfile$interval[1,]), c(3.35333436333173,3.97998278397355))
})

test_that("bmd function computes BMD (extra) correctly for ryegrass model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "inv", display = FALSE)
  resultProfile <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.46370565552042)
  expect_equal(result$bmrScaled[1,1], 7.06180378317457)
  expect_equal(unname(result$interval[1,]), c(1.09762076180767,1.82979054923317))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  # sandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 1.46370565552042)
  expect_equal(resultSandwich$bmrScaled[1,1], 7.06180378317457)
  expect_equal(unname(resultSandwich$interval[1,]), c(1.20293156219483,1.72447974884601))
  # inv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 1.46370565552042)
  expect_equal(resultInv$bmrScaled[1,1], 7.06180378317457)
  expect_equal(unname(resultInv$interval[1,]), c(1.0984899569084,1.78450772463462))
  # profile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 1.46370565552042)
  expect_equal(resultProfile$bmrScaled[1,1], 7.06180378317457)
  expect_equal(unname(resultProfile$interval[1,]), c(1.19113157178988,1.77626232094376))
})

test_that("bmd function computes BMD (relative) correctly for ryegrass model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", interval = "inv", display = FALSE)
  resultProfile <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.49902599632103)
  expect_equal(result$bmrScaled[1,1], 7.01366246432993)
  expect_equal(unname(result$interval[1,]), c(1.169727030614,1.82832496202805))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](0)*0.9))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  # sandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 1.49902599632103)
  expect_equal(resultSandwich$bmrScaled[1,1], 7.01366246432993)
  expect_equal(unname(resultSandwich$interval[1,]), c(1.26495561665339,1.73309637598867))
  # inv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 1.49902599632103)
  expect_equal(resultInv$bmrScaled[1,1], 7.01366246432993)
  expect_equal(unname(resultInv$interval[1,]), c(1.22042713497344,1.85131632264505))
  # profile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 1.49902599632103)
  expect_equal(resultProfile$bmrScaled[1,1], 7.01366246432993)
  expect_equal(unname(resultProfile$interval[1,]), c(1.21827509151193,1.81632858455702))
})

test_that("bmd function computes BMD (added) correctly for ryegrass model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "added", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "added", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "added", backgType = "modelBased", interval = "inv", display = FALSE)
  resultProfile <- bmd(object0, bmr = 0.1, def = "added", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 0.728443033284576)
  expect_equal(result$bmrScaled[1,1], 7.69295829369992)
  expect_equal(unname(result$interval[1,]), c(0.430465534681424,1.02642053188773))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](0)-0.1))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  # sandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 0.728443033284576)
  expect_equal(resultSandwich$bmrScaled[1,1], 7.69295829369992)
  expect_equal(unname(resultSandwich$interval[1,]), c(0.488754570273184,0.968131496295969))
  # inv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 0.728443033284576)
  expect_equal(resultInv$bmrScaled[1,1], 7.69295829369992)
  expect_equal(unname(resultInv$interval[1,]), c(0.538085589755732,1.22374370219594))
  # profile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 0.728443033284576)
  expect_equal(resultProfile$bmrScaled[1,1], 7.69295829369992)
  expect_equal(unname(resultProfile$interval[1,]), c(0.491324061645358,1.03306247242863))
})

test_that("bmd function computes BMD (hybridAdd with hybridSD background) correctly for ryegrass model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, interval = "inv", display = FALSE)
  resultProfile <- bmd(object0, bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, interval = "sandwich", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.21255236145362)
  expect_equal(result$bmrScaled[1,1], 7.35717345395457)
  expect_equal(unname(result$interval[1,]), c(0.873582617250351,1.55152210565689))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  # sandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 1.21255236145362)
  expect_equal(resultSandwich$bmrScaled[1,1], 7.35717345395457)
  expect_equal(unname(resultSandwich$interval[1,]), c(0.963502821126824,1.46160190178042))
  # inv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 1.21255236145362)
  expect_equal(resultInv$bmrScaled[1,1], 7.35717345395457)
  expect_equal(unname(resultInv$interval[1,]), c(0.840460235800311,1.67657493370587))
  # profile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 1.21255236145362)
  expect_equal(resultProfile$bmrScaled[1,1], 7.35717345395457)
  expect_equal(unname(resultProfile$interval[1,]), c(0.963502821126824,1.46160190178042))
})

test_that("bmd function computes BMD (hybridExc with hybridSD background) correctly for ryegrass model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, interval = "inv", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.20672107998472)
  expect_equal(result$bmrScaled[1,1], 7.36302788484143)
  expect_equal(unname(result$interval[1,]), c(0.867947011182731,1.5454951487867))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  # sandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 1.20672107998472)
  expect_equal(resultSandwich$bmrScaled[1,1], 7.36302788484143)
  expect_equal(unname(resultSandwich$interval[1,]), c(0.95757586537566,1.45586629459377))
  # inv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 1.20672107998472)
  expect_equal(resultInv$bmrScaled[1,1], 7.36302788484143)
  expect_equal(unname(resultInv$interval[1,]), c(0.829285739324717,1.6744386486465))
})

test_that("bmd function computes BMD (hybridExc with hybridPercentile background) correctly for ryegrass model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, sandwich.vcov = TRUE, display = FALSE)
  # resultInv <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, interval = "inv", display = FALSE)
  # resultProfile <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, interval = "profile", profileGridSize = 10, display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.06888690340628)
  expect_equal(result$bmrScaled[1,1], 7.4880773398374)
  expect_equal(unname(result$interval[1,]), c(0.73651603354879,1.40125777326377))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  # sandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 1.06888690340628)
  expect_equal(resultSandwich$bmrScaled[1,1], 7.4880773398374)
  expect_equal(unname(resultSandwich$interval[1,]), c(0.818465022872052,1.31930878394051))
})

test_that("bmd function computes BMD (hybridExc with absolute background) correctly for ryegrass model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "absolute", backg = 7, display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "absolute", backg = 7, sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "absolute", backg = 7, interval = "inv", display = FALSE)
  # resultProfile <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, interval = "profile", profileGridSize = 10, display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.02463261154772)
  expect_equal(result$bmrScaled[1,1], 7.52286309492537)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(0.742468238743962,1.30679698435148))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 1.02463261154772)
  expect_equal(resultSandwich$bmrScaled[1,1], 7.52286309492537)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(0.779180809291058,1.27008441380438))
  
  # resultInv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 1.02463261154772)
  expect_equal(resultInv$bmrScaled[1,1], 7.52286309492537)
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_equal(unname(resultInv$interval[1,]), c(0,1.52216300113826))
  
})

test_that("bmd function computes BMD (relative) with log-transformed response correctly for ryegrass model", {
  object0 <- drm(log(rootl) ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", sandwich.vcov = TRUE, display = FALSE)
  expect_error(bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", interval = "inv", display = FALSE),
               "inverse regression interval not available for transformed response.")
  resultProfile <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 0.804218529940602)
  expect_equal(result$bmrScaled[1,1], 2.00113913603417)
  expect_equal(unname(result$interval[1,]), c(0.334333144874303,1.2741039150069))
  expect_equal(result$bmrScaled[1], log(exp(drop(object0$curve[[1]](0)))*0.9))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  # sandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 0.804218529940602)
  expect_equal(resultSandwich$bmrScaled[1,1], 2.00113913603417)
  expect_equal(unname(resultSandwich$interval[1,]), c(0.373462062730365,1.23497499715084))
  # profile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 0.804218529940602)
  expect_equal(resultProfile$bmrScaled[1,1], 2.00113913603417)
  expect_equal(unname(resultProfile$interval[1,]), c(0.44414568898271,1.28705989275751))
})

test_that("bmd function computes BMD (relative) with square root-transformed response correctly for ryegrass model", {
  object0 <- drm(sqrt(rootl) ~ conc, data = drcData::ryegrass, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", sandwich.vcov = TRUE, display = FALSE)
  expect_error(bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", interval = "inv", display = FALSE),
               "inverse regression interval not available for transformed response.")
  resultProfile <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.29590294092622)
  expect_equal(result$bmrScaled[1,1], 2.66331630265024)
  expect_equal(unname(result$interval[1,]), c(0.905813486374546,1.68599239547789))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  # sandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 1.29590294092622)
  expect_equal(resultSandwich$bmrScaled[1,1], 2.66331630265024)
  expect_equal(unname(resultSandwich$interval[1,]), c(1.0388236576693,1.55298222418313))
  # profile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 1.29590294092622)
  expect_equal(resultProfile$bmrScaled[1,1], 2.66331630265024)
  expect_equal(unname(resultProfile$interval[1,]), c(0.970730080509468,1.67207902698738))
})


test_that("bmd function output remains consistent", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  
  snapshot_data <- list(
    Results = as.list(result$Results),
    bmrScaled = as.list(result$bmrScaled),
    interval = as.list(result$interval),
    SE = as.list(result$SE)
  )
  
  # Store a snapshot of the entire result object
  expect_snapshot_value(snapshot_data, style = "deparse")
})



# Ryegrass hormesis -------------------------------------------------------

test_that("bmd function computes BMD (point) correctly for ryegrass hormesis model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmd(object0, bmr = 3.2, def = "point", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 3.63196270261622)
  expect_equal(result$bmrScaled[1,1], 3.2)
  expect_equal(unname(result$interval[1,]), c(3.26967714112966,3.99424826410278))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])), tolerance = 1e-4)
})

test_that("bmd function computes BMD (extra) correctly for ryegrass hormesis model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.53181216366916)
  expect_equal(result$bmrScaled[1,1], 7.00907118396823)
  expect_equal(unname(result$interval[1,]), c(1.14799694948346,1.91562737785485))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])), tolerance = 1e-4)
})

test_that("bmd function computes BMD (relative) correctly for ryegrass hormesis model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.55704870290614)
  expect_equal(result$bmrScaled[1,1], 6.96755641076537)
  expect_equal(unname(result$interval[1,]), c(1.21211587595366,1.90198152985862))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](0)*0.9))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])), tolerance = 1e-4)
})

test_that("bmd function computes BMD (added) correctly for ryegrass hormesis model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmd(object0, bmr = 0.1, def = "added", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.02789452211147)
  expect_equal(result$bmrScaled[1,1], 7.64172934529486)
  expect_equal(unname(result$interval[1,]), c(0.31475182419264,1.7410372200303))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](0)-0.1))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])), tolerance = 1e-4)
})

test_that("bmd function computes BMD (hybridAdd with hybridSD background) correctly for ryegrass hormesis model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.33874261053728)
  expect_equal(result$bmrScaled[1,1], 7.29834344888762)
  expect_equal(unname(result$interval[1,]), c(0.903981308119492,1.77350391295507))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])), tolerance = 1e-4)
})

test_that("bmd function computes BMD (hybridExc with hybridSD background) correctly for ryegrass hormesis model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.33434665932782)
  expect_equal(result$bmrScaled[1,1], 7.30429999406768)
  expect_equal(unname(result$interval[1,]), c(0.897371573938349,1.77132174471728))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])), tolerance = 1e-4)
})

test_that("bmd function computes BMD (hybridExc with hybridPercentile background) correctly for ryegrass hormesis model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.23434488297234)
  expect_equal(result$bmrScaled[1,1], 7.43153058963181)
  expect_equal(unname(result$interval[1,]), c(0.736677604412129,1.73201216153255))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])), tolerance = 1e-4)
})

test_that("bmd function computes BMD (relative) with log-transformed response correctly for ryegrass hormesis model", {
  object0 <- drm(log(rootl) ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.44209342953523)
  expect_equal(result$bmrScaled[1,1], 1.93422445377383)
  expect_equal(unname(result$interval[1,]), c(0.826612118446638,2.05757474062382))
  expect_equal(result$bmrScaled[1], log(exp(drop(object0$curve[[1]](0)))*0.9))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])), tolerance = 1e-4)
})

test_that("bmd function computes BMD (relative) with square root-transformed response correctly for ryegrass hormesis model", {
  object0 <- drm(sqrt(rootl) ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.52253767591639)
  expect_equal(result$bmrScaled[1,1], 2.63448978657361)
  expect_equal(unname(result$interval[1,]), c(1.10162097905078,1.943454372782))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])), tolerance = 1e-4)
})

test_that("bmd function computes correct delta confidence interval for a ryegrass hormesis model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = BC.5())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.53181216366916)
  expect_equal(result$bmrScaled[1,1], 7.00907118396823)
  expect_equal(unname(result$interval[1,]), c(1.14799694948346,1.91562737785485))
})

test_that("bmd function computes correct profile confidence interval for a ryegrass hormesis model", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = BC.5())

  invisible(
    capture.output(
      result <- suppressWarnings(bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "profile",
                                     profileGridSize = 10, display = FALSE, profileProgressInfo = FALSE))))


  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 1.53181216366916)
  expect_equal(result$bmrScaled[1,1], 7.00907118396823)
  expect_equal(unname(result$interval[1,]), c(1.36199561066341,1.71658644037318))
})


# TCDD results ------------------------------------------------------------

test_that("bmd function computes BMD (point) correctly for TCDD model", {
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- bmd(object0, bmr = 0.22, def = "point", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.22, def = "point", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  expect_error(bmd(object0, bmr = 0.22, def = "point", backgType = "modelBased", interval = "inv", display = FALSE),
               "Inverse regression not possible for def=point")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 7.77184985530323)
  expect_equal(result$bmrScaled[1,1], 0.22)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(-9.09858347336425,24.6422831839707), tolerance = 1e-4)
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 7.77184985530323)
  expect_equal(resultSandwich$bmrScaled[1,1], 0.22)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(7.69908281451778,7.84461689608868), tolerance = 1e-4)
  
})

test_that("bmd function computes BMD (excess) correctly for TCDD model", {
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- bmd(object0, bmr = 0.05, def = "excess", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.05, def = "excess", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.05, def = "excess", backgType = "modelBased", interval = "inv", display = FALSE)
  # resultProfile <- bmd(object0, bmr = 0.05, def = "excess", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 5.56116921034511)
  expect_equal(result$bmrScaled[1,1], 0.0709522577318265)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(-30.6802498813352,41.8025883020255), tolerance = 1e-4)
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 5.56116921034511)
  expect_equal(resultSandwich$bmrScaled[1,1], 0.0709522577318265)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(5.40992038310958,5.71241803758064), tolerance = 1e-4)
  
  # resultInv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 5.56116921034511)
  expect_equal(resultInv$bmrScaled[1,1], 0.0709522577318265)
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_equal(unname(resultInv$interval[1,]), c(2.58555484830486,14.4423466517269), tolerance = 1e-4)
  
})

test_that("bmd function computes BMD (additional) correctly for TCDD model", {
  object0 <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  
  result <- bmd(object0, bmr = 0.1, def = "additional", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "additional", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "additional", backgType = "modelBased", interval = "inv", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](0)+0.1))
  
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 6.36475841679501)
  expect_equal(result$bmrScaled[1,1], 0.122055008138765)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(-12.8633221463079,25.5928389798979), tolerance = 1e-4)
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 6.36475841679501)
  expect_equal(resultSandwich$bmrScaled[1,1], 0.122055008138765)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(6.28677997505832,6.4427368585317), tolerance = 1e-4)
  
  # resultInv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 6.36475841679501)
  expect_equal(resultInv$bmrScaled[1,1], 0.122055008138765)
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_equal(unname(resultInv$interval[1,]), c(2.91217704808626,15.3343789956633), tolerance = 1e-4)
  
})



# Chlorac results ---------------------------------------------------------

test_that("bmd function computes BMD (point) correctly for Chlorac model", {
  object0 <- drm(num.dead/total ~ conc, weights = total, data = drcData::chlorac, fct = LN.3u(), type = "binomial")
  
  result <- bmd(object0, bmr = 0.22, def = "point", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.22, def = "point", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultProfile <- suppressWarnings(bmd(object0, bmr = 0.22, def = "point", backgType = "modelBased", interval = "profile", display = FALSE))
  resultProfileGrid <- suppressWarnings(bmd(object0, bmr = 0.22, def = "point", backgType = "modelBased", interval = "profileGrid", profileGridSize = 20, profileProgressInfo = FALSE, display = FALSE))
  expect_error(bmd(object0, bmr = 0.22, def = "point", backgType = "modelBased", interval = "inv", display = FALSE),
               "Inverse regression not possible for def=point")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 22.2271415904074)
  expect_equal(result$bmrScaled[1,1], 0.22)
  expect_equal(unname(result$interval[1,]), c(18.292131969707,26.1621512111078))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 22.2271415904074)
  expect_equal(resultSandwich$bmrScaled[1,1], 0.22)
  expect_equal(unname(resultSandwich$interval[1,]), c(21.7141651046908,22.740118076124))
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  
  # resultProfile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 22.2271415904074)
  expect_equal(resultProfile$bmrScaled[1,1], 0.22)
  expect_equal(unname(resultProfile$interval[1,]), c(19.0428644039017,39.9999389661972))
  
  # resultProfileGrid
  expect_true(!is.na(resultProfileGrid$Results[1, "BMD"]))
  expect_equal(resultProfileGrid$Results[1, "BMD"], 22.2271415904074)
  expect_equal(resultProfileGrid$bmrScaled[1,1], 0.22)
  expect_equal(unname(resultProfileGrid$interval[1,]), c(19.2099579670753,28.9585965214678))
})

test_that("bmd function computes BMD (excess) correctly for Chlorac model", {
  object0 <- drm(num.dead/total ~ conc, weights = total, data = drcData::chlorac, fct = LN.3u(), type = "binomial")
  
  result <- bmd(object0, bmr = 0.05, def = "excess", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.05, def = "excess", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.05, def = "excess", backgType = "modelBased", interval = "inv", display = FALSE)
  resultProfile <- bmd(object0, bmr = 0.05, def = "excess", backgType = "modelBased", interval = "profile", display = FALSE)
  resultProfileGrid <- bmd(object0, bmr = 0.05, def = "excess", backgType = "modelBased", interval = "profileGrid", profileGridSize = 20, profileProgressInfo = FALSE, display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 19.7922909180127)
  expect_equal(result$bmrScaled[1,1], 0.144988275510254)
  expect_equal(unname(result$interval[1,]), c(15.1507479878497,24.4338338481757))
  
  # resultSandwich
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 19.7922909180127)
  expect_equal(resultSandwich$bmrScaled[1,1], 0.144988275510254)
  expect_equal(unname(resultSandwich$interval[1,]), c(18.5884627205015,20.9961191155239))
  
  # resultInv
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 19.7922909180127)
  expect_equal(resultInv$bmrScaled[1,1], 0.144988275510254)
  expect_equal(unname(resultInv$interval[1,]), c(17.2530174112536,25.3502935604775))
  
  # resultProfile
  expect_equal(resultProfile$bmrScaled[1,1], drop(object0$curve[[1]](resultProfile$Results[1, "BMD"])))
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 19.7922909180127)
  expect_equal(resultProfile$bmrScaled[1,1], 0.144988275510254)
  expect_equal(unname(resultProfile$interval[1,]), c(16.1933316166892,39.9999273156964))
  
  # resultProfileGrid
  expect_equal(resultProfileGrid$bmrScaled[1,1], drop(object0$curve[[1]](resultProfileGrid$Results[1, "BMD"])))
  expect_true(!is.na(resultProfileGrid$Results[1, "BMD"]))
  expect_equal(resultProfileGrid$Results[1, "BMD"], 19.7922909180127)
  expect_equal(resultProfileGrid$bmrScaled[1,1], 0.144988275510254)
  expect_equal(unname(resultProfileGrid$interval[1,]), c(16.3460864786873,26.9565322744268))
})

test_that("bmd function computes BMD (additional) correctly for Chlorac model", {
  object0 <- drm(num.dead/total ~ conc, weights = total, data = drcData::chlorac, fct = LN.3u(), type = "binomial")
  
  result <- bmd(object0, bmr = 0.1, def = "additional", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "additional", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "additional", backgType = "modelBased", interval = "inv", display = FALSE)
  resultProfile <- suppressWarnings(bmd(object0, bmr = 0.1, def = "additional", backgType = "modelBased", interval = "profile", display = FALSE))
  resultProfileGrid <- suppressWarnings(bmd(object0, bmr = 0.1, def = "additional", backgType = "modelBased", interval = "profileGrid", profileGridSize = 20, profileProgressInfo = FALSE, display = FALSE))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](0)+0.1))
  
  # result
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 21.7026737107011)
  expect_equal(result$bmrScaled[1,1], 0.199987658431846)
  expect_equal(unname(result$interval[1,]), c(17.173324815284,26.2320226061182))
  
  # resultSandwich
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 21.7026737107011)
  expect_equal(resultSandwich$bmrScaled[1,1], 0.199987658431846)
  expect_equal(unname(resultSandwich$interval[1,]), c(20.4561088161456,22.9492386052567))
  
  # resultInv
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 21.7026737107011)
  expect_equal(resultInv$bmrScaled[1,1], 0.199987658431846)
  expect_equal(unname(resultInv$interval[1,]), c(18.8322967064953,26.4970903630817))
  
  # resultProfile
  expect_equal(resultProfile$bmrScaled[1,1], drop(object0$curve[[1]](resultProfile$Results[1, "BMD"])))
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 21.7026737107011)
  expect_equal(resultProfile$bmrScaled[1,1], 0.199987658431846)
  expect_equal(unname(resultProfile$interval[1,]), c(18.2326808009201,39.9999359543284))
  
  # resultProfileGrid
  expect_equal(resultProfileGrid$bmrScaled[1,1], drop(object0$curve[[1]](resultProfileGrid$Results[1, "BMD"])))
  expect_true(!is.na(resultProfileGrid$Results[1, "BMD"]))
  expect_equal(resultProfileGrid$Results[1, "BMD"], 21.7026737107011)
  expect_equal(resultProfileGrid$bmrScaled[1,1], 0.199987658431846)
  expect_equal(unname(resultProfileGrid$interval[1,]), c(18.2876046579125,28.7023777279169))
  
})


# lemna results -----------------------------------------------------------

test_that("bmd function computes BMD (point) correctly for lemna model", {
  object0 <- drm(frond.num ~ conc, data = drcData::lemna, fct = LL.3(), type = "Poisson")
  
  result <- bmd(object0, bmr = 52, def = "point", backgType = "modelBased", display = FALSE)
  # resultProfile <- bmd(object0, bmr = 52, def = "point", backgType = "modelBased", interval = "profile", display = FALSE)
  
    # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 4.35865965537475)
  expect_equal(result$bmrScaled[1,1], 52)
  expect_equal(unname(result$interval[1,]), c(2.43392036088582,6.28339894986368), tolerance = 1e-4)
  # profile
  # expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  # expect_equal(resultProfile$Results[1, "BMD"], 4.53872215823332)
  # expect_equal(resultProfile$bmrScaled[1,1], 52)
  # expect_equal(unname(resultProfile$interval[1,]), c(3.67365836642196,5.5477249041543))
})

test_that("bmd function computes BMD (extra) correctly for lemna model", {
  object0 <- drm(frond.num ~ conc, data = drcData::lemna, fct = LL.3(), type = "Poisson")
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "inv", display = FALSE)
  # resultProfile <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 0.644966972651776)
  expect_equal(result$bmrScaled[1,1], 60.1147293067283)
  expect_equal(unname(result$interval[1,]), c(-0.274435904843214,1.56436985014677), tolerance = 1e-4)
  # inv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 0.644966972651776)
  expect_equal(resultInv$bmrScaled[1,1], 60.1147293067283)
  expect_equal(unname(resultInv$interval[1,]), c(0.224614484009512,1.33997717698247), tolerance = 1e-4)
  # profile
  # expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  # expect_equal(resultProfile$Results[1, "BMD"], 0.752514405241496)
  # expect_equal(resultProfile$bmrScaled[1,1], 59.8605161777106)
  # expect_equal(unname(resultProfile$interval[1,]), c(0.433615751336297,1.25109127057097))
})

test_that("bmd function computes BMD (relative) correctly for lemna model", {
  object0 <- drm(frond.num ~ conc, data = drcData::lemna, fct = LL.3(), type = "Poisson")
  
  result <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", interval = "inv", display = FALSE)
  # resultProfile <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 0.644966972651776)
  expect_equal(result$bmrScaled[1,1], 60.1147293067283)
  expect_equal(unname(result$interval[1,]), c(-0.124767225483307,1.41470117078686))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](0)*0.9))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])), tolerance = 1e-4)
  # inv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 0.644966972651776)
  expect_equal(resultInv$bmrScaled[1,1], 60.1147293067283)
  expect_equal(unname(resultInv$interval[1,]), c(0.206536602960556,2.00749707239478), tolerance = 1e-4)
  # profile
  # expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  # expect_equal(resultProfile$Results[1, "BMD"], 0.752514405241496)
  # expect_equal(resultProfile$bmrScaled[1,1], 59.8605161777106)
  # expect_equal(unname(resultProfile$interval[1,]), c(0.433615751336297,1.25109127057097))
})




# S.alba results ----------------------------------------------------------

test_that("bmd function computes BMD (point) correctly for S.alba model", {
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  
  result <- bmd(object0, bmr = 3.2, def = "point", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(as.numeric(result$Results[, "BMD"]), c(39.4912945056265, 22.1766859356908))
  expect_equal(as.numeric(result$bmrScaled), c(3.2, 3.2))
  expect_equal(as.numeric(result$bmrScaled), diag(drop(object0$curve[[1]](result$Results[, "BMD"]))))
})

test_that("bmd function computes BMD (relative) correctly for S.alba model", {
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.08, def = "relative", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(as.numeric(result$Results[, "BMD"]), c(28.0790872125237, 18.9735396170819))
  expect_equal(as.numeric(result$bmrScaled), drop(object0$curve[[1]](0))*0.92)
  expect_equal(as.numeric(result$bmrScaled), diag(drop(object0$curve[[1]](result$Results[, "BMD"]))))
})

test_that("bmd function computes BMD (hybridExc with hybridSD background) correctly for S.alba model", {
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(as.numeric(result$Results[, "BMD"]), c(28.0253530227688, 19.0291591246355))
  expect_equal(as.numeric(result$bmrScaled), c(3.56714031355548, 3.49717210738578))
  expect_equal(as.numeric(result$bmrScaled), diag(drop(object0$curve[[1]](result$Results[, "BMD"]))))
})


test_that("bmd function computes correct delta confidence interval for a S.alba model", {
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$interval[, "Lower"])))
  expect_true(all(!is.na(result$interval[, "Upper"])))
  expect_equal(result$Results[, "BMDL"], result$interval[, "Lower"])
  expect_equal(as.numeric(result$interval[, "Lower"]), c(15.4192368012224, 14.5843698054813))
  expect_equal(as.numeric(result$interval[, "Upper"]), c(39.8670300945973, 23.5741646123403))
})

test_that("bmd function computes BMD (point) correctly for S.alba model with independently fitted curves", {
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4(), separate = TRUE)
  
  result <- bmd(object0, bmr = 3.2, def = "point", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(as.numeric(result$Results[, "BMD"]), c(39.4771036152176, 22.1659452714774))
  expect_equal(as.numeric(result$bmrScaled), c(3.2, 3.2))
  expect_equal(as.numeric(result$bmrScaled), diag(drop(object0$curve[[1]](result$Results[, "BMD"]))))
})

test_that("bmd function computes BMD (relative) correctly for S.alba model with independently fitted curves", {
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4(), separate = TRUE)
  
  result <- bmd(object0, bmr = 0.08, def = "relative", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(as.numeric(result$Results[, "BMD"]), c(28.0694995307095, 18.9558015915458))
  expect_equal(as.numeric(result$bmrScaled), drop(object0$curve[[1]](0))*(1-0.08))
  expect_equal(as.numeric(result$bmrScaled), diag(drop(object0$curve[[1]](result$Results[, "BMD"]))))
})


test_that("bmd function output remains consistent with model with multiple curves", {
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  
  snapshot_data <- list(
    Results = as.list(result$Results),
    bmrScaled = as.list(result$bmrScaled),
    interval = as.list(result$interval),
    SE = as.list(result$SE)
  )
  
  # Store a snapshot of the entire result object
  expect_snapshot_value(snapshot_data, style = "deparse")
})


# Increasing continuous model results -------------------------------------

test_that("bmd function computes BMD (point) correctly for increasing continuous model", {
  data0 <- data.frame(dose = c(0,0,0,0,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,
                               25,25,25,25,50,50,50,50,100,100,100,100),
                      resp = c(3.3735,4.1836,3.1644,5.5953,4.6321,3.4822,4.7901,5.041,
                               8.4779,7.5968,9.4139,8.292,15.2536,13.6601,16.9997,15.8299,
                               23.2592,24.2193,24.0967,23.8693,28.7189,28.5821,27.8745,25.8106))
  object0 <- drm(resp ~ dose, data = data0, fct = W1.4())
  
  result <- bmd(object0, bmr = 8, def = "point", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 8, def = "point", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultProfile <- bmd(object0, bmr = 8, def = "point", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # delta
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 12.2806139561059)
  expect_equal(result$bmrScaled[1,1], 8)
  expect_equal(unname(result$interval[1,]), c(11.1666678007201,13.3945601114917))
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  # sandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 12.2806139561059)
  expect_equal(resultSandwich$bmrScaled[1,1], 8)
  expect_equal(unname(resultSandwich$interval[1,]), c(11.4480146049373,13.1132133072745))
  # profile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 12.2806139561059)
  expect_equal(resultProfile$bmrScaled[1,1], 8)
  expect_equal(unname(resultProfile$interval[1,]), c(11.3001432356128,13.2937575997779))
})

test_that("bmd function computes BMD (extra) correctly for increasing continuous model", {
  data0 <- data.frame(dose = c(0,0,0,0,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,
                               25,25,25,25,50,50,50,50,100,100,100,100),
                      resp = c(3.3735,4.1836,3.1644,5.5953,4.6321,3.4822,4.7901,5.041,
                               8.4779,7.5968,9.4139,8.292,15.2536,13.6601,16.9997,15.8299,
                               23.2592,24.2193,24.0967,23.8693,28.7189,28.5821,27.8745,25.8106))
  object0 <- drm(resp ~ dose, data = data0, fct = W1.4())
  
  result <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "inv", display = FALSE)
  resultProfile <- bmd(object0, bmr = 0.1, def = "extra", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 10.863484507375)
  expect_equal(result$bmrScaled[1,1], 6.98897815306841)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(9.51445066757787,12.212518347172))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 10.863484507375)
  expect_equal(resultSandwich$bmrScaled[1,1], 6.98897815306841)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(9.71924181563165,12.0077271991183))
  
  # resultInv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 10.863484507375)
  expect_equal(resultInv$bmrScaled[1,1], 6.98897815306841)
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_equal(unname(resultInv$interval[1,]), c(9.56235201149397,12.2278472729667))
  
  # resultProfile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 10.863484507375)
  expect_equal(resultProfile$bmrScaled[1,1], 6.98897815306841)
  expect_equal(resultProfile$bmrScaled[1,1], drop(object0$curve[[1]](resultProfile$Results[1, "BMD"])))
  expect_equal(unname(resultProfile$interval[1,]), c(9.67873308726007,12.0916306755198))
})

test_that("bmd function computes BMD (relative) correctly for increasing continuous model", {
  data0 <- data.frame(dose = c(0,0,0,0,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,
                               25,25,25,25,50,50,50,50,100,100,100,100),
                      resp = c(3.3735,4.1836,3.1644,5.5953,4.6321,3.4822,4.7901,5.041,
                               8.4779,7.5968,9.4139,8.292,15.2536,13.6601,16.9997,15.8299,
                               23.2592,24.2193,24.0967,23.8693,28.7189,28.5821,27.8745,25.8106))
  object0 <- drm(resp ~ dose, data = data0, fct = W1.4())
  
  result <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", interval = "inv", display = FALSE)
  resultProfile <- suppressWarnings(bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", interval = "profile", display = FALSE))
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 6.41534690365137)
  expect_equal(result$bmrScaled[1,1], 4.55693251645292)
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](0)*1.1))
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(4.75431370333604,8.0763801039667))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 6.41534690365137)
  expect_equal(resultSandwich$bmrScaled[1,1], 4.55693251645292)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(5.0856776962869,7.74501611101584))
  
  # resultInv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 6.41534690365137)
  expect_equal(resultInv$bmrScaled[1,1], 4.55693251645292)
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_equal(unname(resultInv$interval[1,]), c(5.33406644107278,8.89693260160958))
  
  # resultProfile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 6.41534690365137)
  expect_equal(resultProfile$bmrScaled[1,1], 4.55693251645292)
  expect_equal(resultProfile$bmrScaled[1,1], drop(object0$curve[[1]](resultProfile$Results[1, "BMD"])))
  expect_equal(unname(resultProfile$interval[1,]), c(4.95263882310625,7.90896048331597))
  
})

test_that("bmd function computes BMD (added) correctly for increasing continuous model", {
  data0 <- data.frame(dose = c(0,0,0,0,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,
                               25,25,25,25,50,50,50,50,100,100,100,100),
                      resp = c(3.3735,4.1836,3.1644,5.5953,4.6321,3.4822,4.7901,5.041,
                               8.4779,7.5968,9.4139,8.292,15.2536,13.6601,16.9997,15.8299,
                               23.2592,24.2193,24.0967,23.8693,28.7189,28.5821,27.8745,25.8106))
  object0 <- drm(resp ~ dose, data = data0, fct = W1.4())
  
  result <- bmd(object0, bmr = 0.1, def = "added", backgType = "modelBased", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "added", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "added", backgType = "modelBased", interval = "inv", display = FALSE)
  resultProfile <- bmd(object0, bmr = 0.1, def = "added", backgType = "modelBased", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_equal(result$bmrScaled[1], drop(object0$curve[[1]](0)+0.1))
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 4.991720782587)
  expect_equal(result$bmrScaled[1,1], 4.24266592404811)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(3.59642761002814,6.38701395514585))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 4.991720782587)
  expect_equal(resultSandwich$bmrScaled[1,1], 4.24266592404811)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(3.89836064703559,6.0850809181384))
  
  # resultInv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 4.991720782587)
  expect_equal(resultInv$bmrScaled[1,1], 4.24266592404811)
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_equal(unname(resultInv$interval[1,]), c(4.19888523653057,7.381384925039))
  
  # resultProfile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 4.991720782587)
  expect_equal(resultProfile$bmrScaled[1,1], 4.24266592404811)
  expect_equal(resultProfile$bmrScaled[1,1], drop(object0$curve[[1]](resultProfile$Results[1, "BMD"])))
  expect_equal(unname(resultProfile$interval[1,]), c(3.78860392475851,6.26909472927828))
})

test_that("bmd function computes BMD (hybridAdd with hybridSD background) correctly for increasing continuous model", {
  data0 <- data.frame(dose = c(0,0,0,0,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,
                               25,25,25,25,50,50,50,50,100,100,100,100),
                      resp = c(3.3735,4.1836,3.1644,5.5953,4.6321,3.4822,4.7901,5.041,
                               8.4779,7.5968,9.4139,8.292,15.2536,13.6601,16.9997,15.8299,
                               23.2592,24.2193,24.0967,23.8693,28.7189,28.5821,27.8745,25.8106))
  object0 <- drm(resp ~ dose, data = data0, fct = W1.4())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, interval = "inv", display = FALSE)
  resultProfile <- bmd(object0, bmr = 0.1, def = "hybridAdd", backgType = "hybridSD", backg = 2, interval = "sandwich", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 7.5727232324071)
  expect_equal(result$bmrScaled[1,1], 5.00849387544795)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(6.01581326247664,9.12963320233756))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 7.5727232324071)
  expect_equal(resultSandwich$bmrScaled[1,1], 5.00849387544795)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(6.34845172462833,8.79699474018587))
  
  # resultInv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 7.5727232324071)
  expect_equal(resultInv$bmrScaled[1,1], 5.00849387544795)
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_equal(unname(resultInv$interval[1,]), c(6.60326087897833,9.815886863626))
  
  # resultProfile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 7.5727232324071)
  expect_equal(resultProfile$bmrScaled[1,1], 5.00849387544795)
  expect_equal(resultProfile$bmrScaled[1,1], drop(object0$curve[[1]](resultProfile$Results[1, "BMD"])))
  expect_equal(unname(resultProfile$interval[1,]), c(6.34845172462833,8.79699474018587))
  
})

test_that("bmd function computes BMD (hybridExc with hybridSD background) correctly for increasing continuous model", {
  data0 <- data.frame(dose = c(0,0,0,0,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,
                               25,25,25,25,50,50,50,50,100,100,100,100),
                      resp = c(3.3735,4.1836,3.1644,5.5953,4.6321,3.4822,4.7901,5.041,
                               8.4779,7.5968,9.4139,8.292,15.2536,13.6601,16.9997,15.8299,
                               23.2592,24.2193,24.0967,23.8693,28.7189,28.5821,27.8745,25.8106))
  object0 <- drm(resp ~ dose, data = data0, fct = W1.4())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, interval = "inv", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 7.54741579253422)
  expect_equal(result$bmrScaled[1,1], 4.99686214933204)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(5.9914922536827,9.10333933138575))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 7.54741579253422)
  expect_equal(resultSandwich$bmrScaled[1,1], 4.99686214933204)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(6.32405959937184,8.7707719856966))
  
  # resultInv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 7.54741579253422)
  expect_equal(resultInv$bmrScaled[1,1], 4.99686214933204)
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_equal(unname(resultInv$interval[1,]), c(6.5795707606519,9.80594198539089))
  
})

test_that("bmd function computes BMD (hybridExc with hybridPercentile background) correctly for increasing continuous model", {
  data0 <- data.frame(dose = c(0,0,0,0,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,
                               25,25,25,25,50,50,50,50,100,100,100,100),
                      resp = c(3.3735,4.1836,3.1644,5.5953,4.6321,3.4822,4.7901,5.041,
                               8.4779,7.5968,9.4139,8.292,15.2536,13.6601,16.9997,15.8299,
                               23.2592,24.2193,24.0967,23.8693,28.7189,28.5821,27.8745,25.8106))
  object0 <- drm(resp ~ dose, data = data0, fct = W1.4())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, sandwich.vcov = TRUE, display = FALSE)
  # resultInv <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, interval = "inv", display = FALSE)
  # resultProfile <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, interval = "profile", profileGridSize = 10, display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 4.54226313905066)
  expect_equal(result$bmrScaled[1,1], 4.19484787609187)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(3.19216176910602,5.8923645089953))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 4.54226313905066)
  expect_equal(resultSandwich$bmrScaled[1,1], 4.19484787609187)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(3.48330453570381,5.60122174239751))
  
})

test_that("bmd function computes BMD (hybridExc with absolute background) correctly for increasing continuous model", {
  data0 <- data.frame(dose = c(0,0,0,0,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,
                               25,25,25,25,50,50,50,50,100,100,100,100),
                      resp = c(3.3735,4.1836,3.1644,5.5953,4.6321,3.4822,4.7901,5.041,
                               8.4779,7.5968,9.4139,8.292,15.2536,13.6601,16.9997,15.8299,
                               23.2592,24.2193,24.0967,23.8693,28.7189,28.5821,27.8745,25.8106))
  object0 <- drm(resp ~ dose, data = data0, fct = W1.4())
  
  result <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "absolute", backg = 6, display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "absolute", backg = 6, sandwich.vcov = TRUE, display = FALSE)
  resultInv <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "absolute", backg = 6, interval = "inv", display = FALSE)
  # resultProfile <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridPercentile", backg = 0.05, interval = "profile", profileGridSize = 10, display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 7.20955143837895)
  expect_equal(result$bmrScaled[1,1], 4.84877512350992)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(6.0288331932397,8.3902696835182))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 7.20955143837895)
  expect_equal(resultSandwich$bmrScaled[1,1], 4.84877512350992)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(6.40498302800912,8.01411984874878))
  
  # resultInv
  expect_true(!is.na(resultInv$Results[1, "BMD"]))
  expect_equal(resultInv$Results[1, "BMD"], 7.20955143837895)
  expect_equal(resultInv$bmrScaled[1,1], 4.84877512350992)
  expect_equal(resultInv$bmrScaled[1,1], drop(object0$curve[[1]](resultInv$Results[1, "BMD"])))
  expect_equal(unname(resultInv$interval[1,]), c(6.23979423183327,9.47523222361585))
  
})

test_that("bmd function computes BMD (relative) with log-transformed response correctly for increasing continuous model", {
  data0 <- data.frame(dose = c(0,0,0,0,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,
                               25,25,25,25,50,50,50,50,100,100,100,100),
                      resp = c(3.3735,4.1836,3.1644,5.5953,4.6321,3.4822,4.7901,5.041,
                               8.4779,7.5968,9.4139,8.292,15.2536,13.6601,16.9997,15.8299,
                               23.2592,24.2193,24.0967,23.8693,28.7189,28.5821,27.8745,25.8106))
  object0 <- drm(log(resp) ~ dose, data = data0, fct = W1.4())
  
  result <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", sandwich.vcov = TRUE, display = FALSE)
  expect_error(bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", interval = "inv", display = FALSE),
               "inverse regression interval not available for transformed response.")
  resultProfile <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "log", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_equal(result$bmrScaled[1], log(exp(drop(object0$curve[[1]](0)))*1.1))
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 5.68205675583797)
  expect_equal(result$bmrScaled[1,1], 1.46639113206374)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(4.26760401810185,7.09650949357409))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 5.68205675583797)
  expect_equal(resultSandwich$bmrScaled[1,1], 1.46639113206374)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(3.99249625714426,7.37161725453168))
  
  # resultProfile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 5.68205675583797)
  expect_equal(resultProfile$bmrScaled[1,1], 1.46639113206374)
  expect_equal(resultProfile$bmrScaled[1,1], drop(object0$curve[[1]](resultProfile$Results[1, "BMD"])))
  expect_equal(unname(resultProfile$interval[1,]), c(4.49150596950358,7.03069153045658))
  
})

test_that("bmd function computes BMD (relative) with square root-transformed response correctly for increasing continuous model", {
  data0 <- data.frame(dose = c(0,0,0,0,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,
                               25,25,25,25,50,50,50,50,100,100,100,100),
                      resp = c(3.3735,4.1836,3.1644,5.5953,4.6321,3.4822,4.7901,5.041,
                               8.4779,7.5968,9.4139,8.292,15.2536,13.6601,16.9997,15.8299,
                               23.2592,24.2193,24.0967,23.8693,28.7189,28.5821,27.8745,25.8106))
  object0 <- drm(sqrt(resp) ~ dose, data = data0, fct = W1.4())
  
  result <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", display = FALSE)
  resultSandwich <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", sandwich.vcov = TRUE, display = FALSE)
  expect_error(bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", interval = "inv", display = FALSE),
               "inverse regression interval not available for transformed response.")
  resultProfile <- bmd(object0, bmr = 0.1, def = "relative", backgType = "modelBased", respTrans = "sqrt", interval = "profile", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], 5.84243115619949)
  expect_equal(result$bmrScaled[1,1], 2.09788179452544)
  expect_equal(result$bmrScaled[1,1], drop(object0$curve[[1]](result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(4.4433007730042,7.24156153939477))
  
  # resultSandwich
  expect_true(!is.na(resultSandwich$Results[1, "BMD"]))
  expect_equal(resultSandwich$Results[1, "BMD"], 5.84243115619949)
  expect_equal(resultSandwich$bmrScaled[1,1], 2.09788179452544)
  expect_equal(resultSandwich$bmrScaled[1,1], drop(object0$curve[[1]](resultSandwich$Results[1, "BMD"])))
  expect_equal(unname(resultSandwich$interval[1,]), c(4.41009310114347,7.2747692112555))
  
  # resultProfile
  expect_true(!is.na(resultProfile$Results[1, "BMD"]))
  expect_equal(resultProfile$Results[1, "BMD"], 5.84243115619949)
  expect_equal(resultProfile$bmrScaled[1,1], 2.09788179452544)
  expect_equal(resultProfile$bmrScaled[1,1], drop(object0$curve[[1]](resultProfile$Results[1, "BMD"])))
  expect_equal(unname(resultProfile$interval[1,]), c(4.64454353972284,7.12686669438367))
  
})





# Decreasing binomial model with multiple curves --------------------------
test_that("bmd function computes BMD (point) correctly for TCDD model", {
  data0 <- data.frame(
    conc = c(0, rep(c(20, 40, 80, 160, 320), 2)),
    total = c(50, rep(20, 5*2)),
    alive = c(47, 14, 11, 6, 9, 6, 19, 11, 8, 5, 3),
    treat = c("Control", rep(c("Treat1", "Treat2"), each = 5))
  )
  
  object0 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                 pmodels = list(~ treat - 1, ~ treat - 1,
                                ~ 1, ~ treat -1),
                 data = data0, fct = W2.4(), type = "binomial",
                 control = drmc(noMessage = TRUE))
  
  result <- bmd(object0, bmr = 0.77, def = "point", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(result$Results[, "BMD"], c(Treat1 = 16.7094280948318, Treat2 = 30.2464214935383))
  expect_equal(result$bmrScaled[,1], c(Treat1 = 0.77, Treat2 = 0.77))
  expect_equal(unname(result$bmrScaled[,1]), diag(object0$curve[[1]](result$Results[, "BMD"])[,2:3]))
  expect_equal(result$interval[, "Lower"], c(Treat1 = 6.3733358375481, Treat2 = 22.2685483133783))
  expect_equal(result$interval[, "Upper"], c(Treat1 = 27.0455203521155, Treat2 = 38.2242946736982))
})

test_that("bmd function computes BMD (excess) correctly for TCDD model", {
  data0 <- data.frame(
    conc = c(0, rep(c(20, 40, 80, 160, 320), 2)),
    total = c(50, rep(20, 5*2)),
    alive = c(47, 14, 11, 6, 9, 6, 19, 11, 8, 5, 3),
    treat = c("Control", rep(c("Treat1", "Treat2"), each = 5))
  )
  
  object0 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                 pmodels = list(~ treat - 1, ~ treat - 1,
                                ~ 1, ~ treat -1),
                 data = data0, fct = W2.4(), type = "binomial",
                 control = drmc(noMessage = TRUE))
  
  result <- bmd(object0, bmr = 0.1, def = "excess", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(result$Results[, "BMD"], c(Treat1 = 12.7945107782873, Treat2 = 25.107888955923))
  expect_equal(result$bmrScaled[,1], c(Treat1 = 0.851586285814034, Treat2 = 0.85158628581403477))
  expect_equal(unname(result$bmrScaled[,1]), diag(object0$curve[[1]](result$Results[, "BMD"])[,2:3]))
  expect_equal(result$interval[, "Lower"], c(Treat1 = 1.21746148052096, Treat2 = 16.6628970797628))
  expect_equal(result$interval[, "Upper"], c(Treat1 = 24.3715600760537, Treat2 = 33.5528808320832))
  
})

test_that("bmd function computes BMD (additional) correctly for TCDD model", {
  data0 <- data.frame(
    conc = c(0, rep(c(20, 40, 80, 160, 320), 2)),
    total = c(50, rep(20, 5*2)),
    alive = c(47, 14, 11, 6, 9, 6, 19, 11, 8, 5, 3),
    treat = c("Control", rep(c("Treat1", "Treat2"), each = 5))
  )
  
  object0 <- drm(alive/total ~ conc, weights = total, curveid = treat, 
                 pmodels = list(~ treat - 1, ~ treat - 1,
                                ~ 1, ~ treat -1),
                 data = data0, fct = W2.4(), type = "binomial",
                 control = drmc(noMessage = TRUE))
  
  result <- bmd(object0, bmr = 0.1, def = "additional", backgType = "modelBased", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # result
  expect_true(all(!is.na(result$Results[, "BMD"])))
  expect_equal(result$Results[, "BMD"], c(Treat1 = 13.0508853789932, Treat2 = 25.4655976138905))
  expect_equal(result$bmrScaled[,1], c(Treat1 = 0.846206984237815, Treat2 = 0.846206984237815))
  expect_equal(unname(result$bmrScaled[,1]), diag(object0$curve[[1]](result$Results[, "BMD"])[,2:3]))
  expect_equal(result$interval[, "Lower"], c(Treat1 = 1.50118005072325, Treat2 = 16.964077339049))
  expect_equal(result$interval[, "Upper"], c(Treat1 = 24.6005907072632, Treat2 = 33.9671178887321))
  
})


# Meta-analytic random effects model --------------------------------------

test_that("bmd function works on drcMMRE object", {
  set.seed(1)
  data0 <- data.frame(x = rep(drcData::ryegrass$conc, 2),
                      y = rep(drcData::ryegrass$rootl, 2) +
                        c(rnorm(n = nrow(drcData::ryegrass), mean = 2, sd = 0.5),
                          rnorm(n = nrow(drcData::ryegrass), mean = 2.7, sd = 0.7)),
                      EXP_ID = rep(as.character(1:2), each = nrow(drcData::ryegrass)))
  
  modMMRE <- drmMMRE(y~x, exp_id = EXP_ID, data = data0, fct = LL.4())
  bmdMMRE <- bmd(modMMRE, bmr = 0.1, backgType = "modelBased", def = "relative", display = FALSE)
  
  expect_true(all(!is.na(bmdMMRE$Results[, "BMD"])))
  expect_equal(bmdMMRE$Results[, "BMD"], 1.66913593445629)
  expect_equal(bmdMMRE$bmrScaled[,1], 9.15712352078559)
  expect_equal(unname(bmdMMRE$bmrScaled[,1]), drop(modMMRE$curve[[1]](bmdMMRE$Results[, "BMD"])))
  expect_equal(bmdMMRE$interval[1,], c(Lower = 1.3166277025622, Upper = 2.02164416635037), tolerance = 1e-4)
  expect_equal(bmdMMRE$SE[,"SE"], 0.214309787885148, tolerance = 1e-4) 
})
