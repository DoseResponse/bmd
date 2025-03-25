# Tests for bmdHetVar function
# - Arguments and structure
#   - Missing arguments
#   - Correct backgType and def accepted
# - Ryegrass model
#   - correct bmd estimate (hybridExc + hybridAdd with hybridPercentile and hybridSD backgType)
# - GiantKelp model
#   - correct bmd estimate (hybridExc + hybridAdd with hybridPercentile and hybridSD backgType)




# Arguments and structure -------------------------------------------------

test_that("bmdHetVar handles missing required arguments", {
  lm_object <- lm(y ~ x, 
                  data = data.frame(x = 0:4, 
                                    y = 1:5 + c(-0.4, 0, 0.2, -0.1, 0.13)))
  ryegrass.W2.4 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = W2.4())
  var.formula0 <- ~ fitted + I(fitted^2)
  ryegrass.W2.4.hetVar <- drmHetVar(ryegrass.W2.4, var.formula0)
  
  expect_error(bmdHetVar(lm_object), 'object must be a dose-response model with a heterogeneous variance structure of class "drcHetVar"')
  expect_error(bmdHetVar(ryegrass.W2.4), 'object must be a dose-response model with a heterogeneous variance structure of class "drcHetVar"')
  expect_error(bmdHetVar(ryegrass.W2.4.hetVar, bmr = 0.1, def = "hybridExc"), 'backgType is missing. Options are "absolute", "hybridSD" or "hybridPercentile"')
  expect_error(bmdHetVar(ryegrass.W2.4.hetVar, bmr = 0.1, def = "hybridExc", backgType = "wrongType"), 'Could not recognize backgType. Options are "absolute", "hybridSD" or "hybridPercentile"')
  expect_error(bmdHetVar(bmr = 0.1, backgType = "hybridPercentile", def = "hybridExc"), 'argument "object" is missing, with no default')
  expect_error(bmdHetVar(ryegrass.W2.4.hetVar, backgType = "hybridPercentile", def = "hybridExc"), 'argument "bmr" needs to be specified as a number between 0 and 1')
  expect_error(bmdHetVar(ryegrass.W2.4.hetVar, bmr = 1.5, backgType = "hybridPercentile", def = "hybridExc"), 'argument "bmr" needs to be specified as a number between 0 and 1')
  expect_error(bmdHetVar(ryegrass.W2.4.hetVar, bmr = 0.1, backgType = "hybridPercentile"), 'def is missing. Options are "hybridExc" or "hybridAdd"')
  expect_error(bmdHetVar(ryegrass.W2.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", def = "wrongDef"), 'Could not recognize def. Options are "hybridExc" or "hybridAdd"')
})





# Ryegrass model ----------------------------------------------------------

test_that("bmdHetVar on Ryegrass model, def = hybridExc, backgtype = hybridPercentile", {
  ryegrass.W2.4 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = W2.4())
  var.formula0 <- ~ fitted + I(fitted^2)
  ryegrass.W2.4.hetVar <- drmHetVar(ryegrass.W2.4, var.formula0)
  set.seed(1)
  result <- bmdHetVar(ryegrass.W2.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c("BMD"=1.27811671431062))
  expect_equal(result$bmrScaled[1,1], 7.54478765441619)
  expect_equal(unname(result$interval[1,]), c(1.12131595768687,1.43679973459246))
})

test_that("bmdHetVar on Ryegrass model, def = hybridExc, backgtype = hybridSD", {
  ryegrass.W2.4 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = W2.4())
  var.formula0 <- ~ fitted + I(fitted^2)
  ryegrass.W2.4.hetVar <- drmHetVar(ryegrass.W2.4, var.formula0)
  set.seed(1)
  result <- bmdHetVar(ryegrass.W2.4.hetVar, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridExc", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c("BMD"=1.45386495506409))
  expect_equal(result$bmrScaled[1,1], 7.30984899101547)
  expect_equal(unname(result$interval[1,]), c(1.26150389977647,1.57147862329125))
})

test_that("bmdHetVar on Ryegrass model, def = hybridAdd, backgtype = hybridPercentile", {
  ryegrass.W2.4 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = W2.4())
  var.formula0 <- ~ fitted + I(fitted^2)
  ryegrass.W2.4.hetVar <- drmHetVar(ryegrass.W2.4, var.formula0)
  set.seed(1)
  result <- bmdHetVar(ryegrass.W2.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c("BMD"=1.2951863924654))
  expect_equal(result$bmrScaled[1,1], 7.52658664968204)
  expect_equal(unname(result$interval[1,]), c(1.13769226756753,1.45126563749746))
})

test_that("bmdHetVar on Ryegrass model, def = hybridAdd, backgtype = hybridSD", {
  ryegrass.W2.4 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = W2.4())
  var.formula0 <- ~ fitted + I(fitted^2)
  ryegrass.W2.4.hetVar <- drmHetVar(ryegrass.W2.4, var.formula0)
  set.seed(1)
  result <- bmdHetVar(ryegrass.W2.4.hetVar, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c("BMD"=1.45814046464533))
  expect_equal(result$bmrScaled[1,1], 7.30288607148892)
  expect_equal(unname(result$interval[1,]), c(1.26409211540589,1.57578126779336))
})



# GiantKelp model ---------------------------------------------------------


test_that("bmdHetVar on GiantKelp model, def = hybridExc, backgtype = hybridPercentile", {
  GiantKelp.LL.4 <- drm(tubeLength ~ dose, data = drcData::GiantKelp, fct = LL.4())
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  GiantKelp.LL.4.hetVar <- drmHetVar(GiantKelp.LL.4, var.formula0)
  set.seed(1)
  result <- bmdHetVar(GiantKelp.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c("BMD"=5.13707358300282))
  expect_equal(result$bmrScaled[1,1], 17.3098274841786)
  expect_equal(unname(result$interval[1,]), c(0.888110687020647,14.4030845590698))
})

test_that("bmdHetVar on GiantKelp model, def = hybridExc, backgtype = hybridSD", {
  GiantKelp.LL.4 <- drm(tubeLength ~ dose, data = drcData::GiantKelp, fct = LL.4())
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  GiantKelp.LL.4.hetVar <- drmHetVar(GiantKelp.LL.4, var.formula0)
  set.seed(1)
  result <- bmdHetVar(GiantKelp.LL.4.hetVar, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridExc", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c("BMD"=10.0346000078297))
  expect_equal(result$bmrScaled[1,1], 16.472304344928)
  expect_equal(unname(result$interval[1,]), c(3.40106451148093,20.503469300298))
})

test_that("bmdHetVar on GiantKelp model, def = hybridAdd, backgtype = hybridPercentile", {
  GiantKelp.LL.4 <- drm(tubeLength ~ dose, data = drcData::GiantKelp, fct = LL.4())
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  GiantKelp.LL.4.hetVar <- drmHetVar(GiantKelp.LL.4, var.formula0)
  set.seed(1)
  result <- bmdHetVar(GiantKelp.LL.4.hetVar, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c("BMD"=5.57094720723732))
  expect_equal(result$bmrScaled[1,1], 17.2358870698386)
  expect_equal(unname(result$interval[1,]), c(1.02451690713505,15.0145292097429))
})

test_that("bmdHetVar on GiantKelp model, def = hybridAdd, backgtype = hybridSD", {
  GiantKelp.LL.4 <- drm(tubeLength ~ dose, data = drcData::GiantKelp, fct = LL.4())
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  GiantKelp.LL.4.hetVar <- drmHetVar(GiantKelp.LL.4, var.formula0)
  set.seed(1)
  result <- bmdHetVar(GiantKelp.LL.4.hetVar, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c("BMD"=10.1611348555932))
  expect_equal(result$bmrScaled[1,1], 16.4508957731788)
  expect_equal(unname(result$interval[1,]), c(3.49771946179717,20.6442510937985))
})



