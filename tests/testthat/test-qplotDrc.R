# Tests for qplotDrc function
# TCDD model
# S.alba model


# TCDD model --------------------------------------------------------------
test_that("qplotDrc returns a ggplot object", {
  object.TCDD <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  p <- qplotDrc(object.TCDD)
  expect_s3_class(p, "ggplot")
})

test_that("qplotDrc handles log transformation", {
  object.TCDD <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  p <- qplotDrc(object.TCDD, xtrans = "log")
  expect_s3_class(p, "ggplot")
})

test_that("qplotDrc handles pseudo_log transformation", {
  object.TCDD <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  p <- qplotDrc(object.TCDD, xtrans = "pseudo_log")
  expect_s3_class(p, "ggplot")
})



test_that("qplotDrc fails when level is out of range", {
  object.TCDD <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  expect_error(qplotDrc(object.TCDD, level = "invalid_level"), "Nothing to plot")
})

test_that("qplotDrc: add argument works", {
  object.TCDD <- drm(incidence/total ~ conc, weights = total, fct = LL.4(), data = drcData::TCDD, type = "binomial")
  p <- qplotDrc(object.TCDD) +
    qplotDrc(object.TCDD, type = "confidence", add = TRUE)$confBandLayer
  expect_s3_class(p, "ggplot")
})




# S.alba model ------------------------------------------------------------

test_that("qplotDrc handles different plot types", {
  object.S.alba <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  types <- c("average", "all", "bars", "none", "obs", "confidence")
  for (t in types) {
    p <- qplotDrc(object.S.alba, type = t)
    expect_s3_class(p, "ggplot")
  }
})


test_that("qplotDrc handles colours", {
  object.S.alba <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  p <- qplotDrc(object.S.alba, col = TRUE)
  expect_s3_class(p, "ggplot")
})



