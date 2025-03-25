# Tests for bmdOrdinal function
# - Arguments and structure
#   - Missing arguments
#   - Correct backgType and def accepted
# - guthionS model
#   - correct bmd estimate (all definitions)
#   - delta and bootstrap intervals




# Arguments and structure -------------------------------------------------

test_that("bmdOrdinal handles missing required arguments",{
  guthionS <- subset(drcData::guthion, trt == "S")
  
  expect_error(drmOrdinal(weights = "total", dose = "dose", data = guthionS, fct = LL.2()), 'argument "levels" is missing, with no default')
  expect_error(drmOrdinal(levels = c("alive", "moribund", "dead"), dose = "dose", data = guthionS, fct = LL.2()), 'argument "weights" is missing, with no default')
  expect_error(drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose",fct = LL.2()), 'argument "data" is missing, with no default')
  expect_error(drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS), 'argument "fct" is missing, with no default')
})




# guthionS model ----------------------------------------------------------

test_that("bmdOrdinal computes correct bmd estimates (def = point) for gutionS model", {
  guthionS <- subset(drcData::guthion, trt == "S")
  
  guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = LL.2())
  
  # checks
  expect_equal(guthionS.LL$levels, c("alive", "moribund", "dead"))
  expect_equal(guthionS.LL$levelsMerged, list("moribund+dead", "dead"))
  expect_equal(guthionS.LL$dose, "dose")
  expect_equal(guthionS.LL$weights, "total")
  expect_null(guthionS.LL$blocks)
  
  expect_equal(guthionS.LL$pFun(0), c("alive" = 1, "moribund" = 0, "dead" = 0))
  expect_equal(guthionS.LL$pFun(10), c("alive" = 0.995187445898751, "moribund" = 0, "dead" = 0.00725091623314945))
  expect_equal(guthionS.LL$pFun(20), c("alive" = 0.909311643501763, "moribund" = 0.000163464702208557, "dead" = 0.0905248917960284))
  expect_equal(guthionS.LL$pFun(30), c("alive" = 0.630625814086567, "moribund" = 0.0548943439739693, "dead" = 0.314479841939464))
  expect_equal(guthionS.LL$pFun(40), c("alive" = 0.327132011789939, "moribund" = 0.0972367191475815, "dead" = 0.57563126906248))
  expect_equal(guthionS.LL$pFun(45), c("alive" = 0.225228290238974, "moribund" = 0.0958704725326432, "dead" = 0.678901237228382))
  
  expect_equal(AIC(guthionS.LL), 222.869282843461)
})
