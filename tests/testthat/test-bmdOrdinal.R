# Tests for drmOrdinal function
# - Arguments and structure
#   - Missing arguments
#   - Correct backgType and def accepted
# - guthionS model
#   - correct bmd estimate (all definitions)
#   - delta and bootstrap intervals




# Arguments and structure -------------------------------------------------

## TO BE ADDED




# guthionS model ----------------------------------------------------------

test_that("bmdOrdinal computes correct bmd estimates (def = point) for gutionS model", {
  guthionS <- subset(drcData::guthion, trt == "S")
  
  guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = LL.2())
  result <- bmdOrdinal(guthionS.LL, bmr = 0.2, backgType = "modelBased", def = "point", display = FALSE)
  
  # checks
  expect_true(all(!is.na(result$Results[,1])))
  expect_equal(result$Results[,1], c("moribund+dead" = 24.6851764586311, "dead" = 25.5366473289933), tolerance = 1e-4)
  expect_equal(result$Results[,2], c("moribund+dead" = 21.6860301841392, "dead" = 22.153171404582), tolerance = 1e-4)
  expect_equal(result$interval[,2], c("moribund+dead" = 27.684322733123, "dead" = 28.9201232534045), tolerance = 1e-4)
})

test_that("bmdOrdinal computes correct bmd estimates (def = excess) for gutionS model", {
  guthionS <- subset(drcData::guthion, trt == "S")
  
  guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = LL.2())
  result <- bmdOrdinal(guthionS.LL, bmr = 0.1, backgType = "modelBased", def = "excess", display = FALSE)
  
  # checks
  expect_true(all(!is.na(result$Results[,1])))
  expect_equal(result$Results[,1], c("moribund+dead" = 20.501035524755, "dead" = 20.5924938340352), tolerance = 1e-4)
  expect_equal(result$Results[,2], c("moribund+dead" = 17.1148425800741, "dead" = 16.6957341161554), tolerance = 1e-4)
  expect_equal(result$interval[,2], c("moribund+dead" = 23.887228469436, "dead" = 24.4892535519151), tolerance = 1e-4)
})

test_that("bmdOrdinal computes correct bmd estimates (def = additional) for gutionS model", {
  guthionS <- subset(drcData::guthion, trt == "S")
  
  guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = LL.2())
  result <- bmdOrdinal(guthionS.LL, bmr = 0.1, backgType = "modelBased", def = "additional", display = FALSE)
  
  # checks
  expect_true(all(!is.na(result$Results[,1])))
  expect_equal(result$Results[,1], c("moribund+dead" = 20.501035524755, "dead" = 20.5924938340352), tolerance = 1e-4)
  expect_equal(result$Results[,2], c("moribund+dead" = 17.1148425800741, "dead" = 16.6957341161554), tolerance = 1e-4)
  expect_equal(result$interval[,2], c("moribund+dead" = 23.887228469436, "dead" = 24.4892535519151), tolerance = 1e-4)
})

test_that("bmdOrdinal computes correct confidence intervals for gutionS model", {
  guthionS <- subset(drcData::guthion, trt == "S")
  
  guthionS.LL <- drmOrdinal(levels = c("alive", "moribund", "dead"), weights = "total", dose = "dose", data = guthionS, fct = LL.2())
  result.delta <- bmdOrdinal(guthionS.LL, bmr = 0.1, backgType = "modelBased", def = "excess", display = FALSE)
  result.sandwich <- bmdOrdinal(guthionS.LL, bmr = 0.1, backgType = "modelBased", def = "excess", interval = "sandwich", display = FALSE)
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result.boot <- bmdOrdinal(guthionS.LL, bmr = 0.1, backgType = "modelBased", def = "excess", interval = "bootstrap", R = 50, display = FALSE, progressInfo = FALSE)
  
  # checks
  expect_true(all(!is.na(result.delta$Results[,1])))
  expect_equal(result.delta$Results[,1], c("moribund+dead" = 20.501035524755, "dead" = 20.5924938340352), tolerance = 1e-4)
  expect_equal(result.sandwich$Results[,1], c("moribund+dead" = 20.501035524755, "dead" = 20.5924938340352), tolerance = 1e-4)
  expect_equal(result.boot$Results[,1], c("moribund+dead" = 20.501035524755, "dead" = 20.5924938340352), tolerance = 1e-4)
  
  expect_equal(result.delta$interval[,1], c("moribund+dead" = 17.1148425800741, "dead" = 16.6957341161554), tolerance = 1e-4)
  expect_equal(result.delta$interval[,2], c("moribund+dead" = 23.887228469436, "dead" = 24.4892535519151), tolerance = 1e-4)
  
  expect_equal(result.sandwich$interval[,1], c("moribund+dead" = 17.712364940837, "dead" = 19.5497322478028), tolerance = 1e-4)
  expect_equal(result.sandwich$interval[,2], c("moribund+dead" = 23.289706108673, "dead" = 21.6352554202677), tolerance = 1e-4)
  
  expect_equal(result.boot$interval[,1], c("moribund+dead" = 17.7581615954546, "dead" = 16.8551004638107), tolerance = 1e-4)
  expect_equal(result.boot$interval[,2], c("moribund+dead" = 24.1801039423246, "dead" = 24.3514247473655), tolerance = 1e-4)
})
