# Tests for qplotBmd function
# TCDD model
# S.alba model


# Arguments and structure -------------------------------------------------

test_that("qplotBmd handles wrong objects and arguments", {
  object0 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  object0.LN <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LN.4())
  resultMA <- bmdMA(list(object0, object0.LN), modelWeights = "AIC", bmr = 3.2, def = "point", backgType = "modelBased", type = "Kang", display = FALSE)
  
  expect_error(qplotBmd(object0), 'qplotBmd only works for plotting objects of type "bmd"')
  expect_error(qplotBmd(resultMA), 'qplotBmd does not for for model-averaged BMD')
})


# TCDD model --------------------------------------------------------------
test_that("qplotBmd returns a ggplot object", {
  object.TCDD <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  bmd.TCDD <- bmd(object.TCDD, bmr = 0.05, def = "excess", backgType = "modelBased", display = FALSE)
  
  p <- qplotBmd(bmd.TCDD)
  expect_s3_class(p, "ggplot")
})


test_that("qplotBmd: add argument works", {
  object.TCDD <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  bmd.TCDD <- bmd(object.TCDD, bmr = 0.05, def = "excess", backgType = "modelBased", sandwich.vcov = TRUE, display = FALSE)
  
  p <- qplotDrc(object.TCDD) +
    qplotDrc(object.TCDD, type = "confidence", add = TRUE)$confBandLayer +
    qplotBmd(bmd.TCDD, add = TRUE)
  expect_s3_class(p, "ggplot")
})




# S.alba model ------------------------------------------------------------

test_that("qplotBmd handles different plot types", {
  object.S.alba <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  bmd.S.alba <- bmd(object.S.alba, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, display = FALSE)
  
  types <- c("average", "all", "bars", "none", "obs", "confidence")
  for (t in types) {
    p <- qplotDrc(object.S.alba, type = t) + qplotBmd(bmd.S.alba, add = TRUE)
    expect_s3_class(p, "ggplot")
  }
})


test_that("qplotBmd handles colours", {
  object.S.alba <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  bmd.S.alba <- bmd(object.S.alba, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, display = FALSE)
  
  p <- qplotBmd(bmd.S.alba, col = TRUE)
  expect_s3_class(p, "ggplot")
})



