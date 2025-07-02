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
  ryegrass.W2.4.hetVar <- drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W2.4())
  
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
  var.formula0 <- ~ fitted + I(fitted^2)
  object0 <- drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W2.4())
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  resultSemiParametric <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 50, level = 0.95, bootType = "semiparametric", progressInfo = FALSE, display = FALSE)
  resultParametric <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 50, level = 0.95, bootType = "parametric", progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c(BMD=1.22216776066126))
  expect_equal(result$bmrScaled[1,1], 7.53945959231123)
  expect_equal(result$bmrScaled[1,1], unname(object0$curve(result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(1.04983687382092,1.41290802464387), tolerance = 1e-1)
  
  # resultSemiParametric
  expect_true(!is.na(resultSemiParametric$Results[1, "BMD"]))
  expect_equal(resultSemiParametric$Results[1, "BMD"], c(BMD=1.22216776066126))
  expect_equal(resultSemiParametric$bmrScaled[1,1], 7.53945959231123)
  expect_equal(resultSemiParametric$bmrScaled[1,1], unname(object0$curve(resultSemiParametric$Results[1, "BMD"])))
  expect_equal(unname(resultSemiParametric$interval[1,]), c(0.950640020523221,1.41523185031235), tolerance = 1e-1)
  
  # resultParametric
  expect_true(!is.na(resultParametric$Results[1, "BMD"]))
  expect_equal(resultParametric$Results[1, "BMD"], c(BMD=1.22216776066126))
  expect_equal(resultParametric$bmrScaled[1,1], 7.53945959231123)
  expect_equal(resultParametric$bmrScaled[1,1], unname(object0$curve(resultParametric$Results[1, "BMD"])))
  expect_equal(unname(resultParametric$interval[1,]), c(1.00050130374069,1.44861779926571), tolerance = 1e-1)
  
})

test_that("bmdHetVar on Ryegrass model, def = hybridExc, backgtype = hybridSD", {
  var.formula0 <- ~ fitted + I(fitted^2)
  object0 <- drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W2.4())
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridExc", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  resultParametric <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridExc", R = 50, level = 0.95, bootType = "parametric", progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c(BMD=1.40198831069816))
  expect_equal(result$bmrScaled[1,1], 7.29990540347001)
  expect_equal(result$bmrScaled[1,1], unname(object0$curve(result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(1.21029877407248,1.5675419352485), tolerance = 1e-1)
  
  # resultParametric
  expect_true(!is.na(resultParametric$Results[1, "BMD"]))
  expect_equal(resultParametric$Results[1, "BMD"], c(BMD=1.40198831069816))
  expect_equal(resultParametric$bmrScaled[1,1], 7.29990540347001)
  expect_equal(resultParametric$bmrScaled[1,1], unname(object0$curve(resultParametric$Results[1, "BMD"])))
  expect_equal(unname(resultParametric$interval[1,]), c(1.17597981534217,1.63561245267409), tolerance = 1e-1)
  
})

test_that("bmdHetVar on Ryegrass model, def = hybridAdd, backgtype = hybridPercentile", {
  var.formula0 <- ~ fitted + I(fitted^2)
  object0 <- drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W2.4())
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  resultParametric <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridAdd", R = 50, level = 0.95, bootType = "parametric", progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c(BMD=1.23980768008069))
  expect_equal(result$bmrScaled[1,1], 7.52040683311043)
  expect_equal(result$bmrScaled[1,1], unname(object0$curve(result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(1.06457123469667,1.42648905128196), tolerance = 1e-1)
  
  # resultParametric
  expect_true(!is.na(resultParametric$Results[1, "BMD"]))
  expect_equal(resultParametric$Results[1, "BMD"], c(BMD=1.23980768008069))
  expect_equal(resultParametric$bmrScaled[1,1], 7.52040683311043)
  expect_equal(resultParametric$bmrScaled[1,1], unname(object0$curve(resultParametric$Results[1, "BMD"])))
  expect_equal(unname(resultParametric$interval[1,]), c(1.03249045631134,1.51113175875138), tolerance = 1e-1)
  
})

test_that("bmdHetVar on Ryegrass model, def = hybridAdd, backgtype = hybridSD", {
  var.formula0 <- ~ fitted + I(fitted^2)
  object0 <- drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W2.4())
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  resultParametric <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridAdd", R = 50, level = 0.95, bootType = "parametric", progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c(BMD=1.40629441578573))
  expect_equal(result$bmrScaled[1,1], 7.29301416210028)
  expect_equal(result$bmrScaled[1,1], unname(object0$curve(result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(1.21501660743054,1.57188381181681), tolerance = 1e-1)
  
  # resultParametric
  expect_true(!is.na(resultParametric$Results[1, "BMD"]))
  expect_equal(resultParametric$Results[1, "BMD"], c(BMD=1.40629441578573))
  expect_equal(resultParametric$bmrScaled[1,1], 7.29301416210028)
  expect_equal(resultParametric$bmrScaled[1,1], unname(object0$curve(resultParametric$Results[1, "BMD"])))
  expect_equal(unname(resultParametric$interval[1,]), c(1.18001446789916,1.63934271919724), tolerance = 1e-1)
  
})



# GiantKelp model ---------------------------------------------------------


test_that("bmdHetVar on GiantKelp model, def = hybridExc, backgtype = hybridPercentile", {
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  object0 <- drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LL.4())
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  resultSemiParametric <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 50, level = 0.95, bootType = "semiparametric", progressInfo = FALSE, display = FALSE)
  resultParametric <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridExc", R = 50, level = 0.95, bootType = "parametric", progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c(BMD=5.62878374597309))
  expect_equal(result$bmrScaled[1,1], 17.4208624811968)
  expect_equal(result$bmrScaled[1,1], unname(object0$curve(result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(1.04968526838944,15.522684535877), tolerance = 1e-1)
  
  # resultSemiParametric
  expect_true(!is.na(resultSemiParametric$Results[1, "BMD"]))
  expect_equal(resultSemiParametric$Results[1, "BMD"], c(BMD=5.62878374597309))
  expect_equal(resultSemiParametric$bmrScaled[1,1], 17.4208624811968)
  expect_equal(resultSemiParametric$bmrScaled[1,1], unname(object0$curve(resultSemiParametric$Results[1, "BMD"])))
  expect_equal(unname(resultSemiParametric$interval[1,]), c(3.83490020339583,10.0947390719072), tolerance = 1e-1)
  
  # resultParametric
  expect_true(!is.na(resultParametric$Results[1, "BMD"]))
  expect_equal(resultParametric$Results[1, "BMD"], c(BMD=5.62878374597309))
  expect_equal(resultParametric$bmrScaled[1,1], 17.4208624811968)
  expect_equal(resultParametric$bmrScaled[1,1], unname(object0$curve(resultParametric$Results[1, "BMD"])))
  expect_equal(unname(resultParametric$interval[1,]), c(3.18753866045103,10.1753755402802), tolerance = 1e-1)
})

test_that("bmdHetVar on GiantKelp model, def = hybridExc, backgtype = hybridSD", {
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  object0 <- drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LL.4())
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridExc", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  resultParametric <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridExc", R = 50, level = 0.95, bootType = "parametric", progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c(BMD=10.3172737799034))
  expect_equal(result$bmrScaled[1,1], 16.5889266156056)
  expect_equal(result$bmrScaled[1,1], unname(object0$curve(result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(3.34558983991111,21.1813243422884), tolerance = 1e-1)
  
  # resultParametric
  expect_true(!is.na(resultParametric$Results[1, "BMD"]))
  expect_equal(resultParametric$Results[1, "BMD"], c(BMD=10.3172737799034))
  expect_equal(resultParametric$bmrScaled[1,1], 16.5889266156056)
  expect_equal(resultParametric$bmrScaled[1,1], unname(object0$curve(resultParametric$Results[1, "BMD"])))
  expect_equal(unname(resultParametric$interval[1,]), c(5.9069620145796,14.6332619325219), tolerance = 1e-1)
  
})

test_that("bmdHetVar on GiantKelp model, def = hybridAdd, backgtype = hybridPercentile", {
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  object0 <- drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LL.4())
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridPercentile", backg = 0.1, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c(BMD=6.05249924528408))
  expect_equal(result$bmrScaled[1,1], 17.3478591564324)
  expect_equal(result$bmrScaled[1,1], unname(object0$curve(result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(1.23069860312066,16.0998877338217), tolerance = 1e-1)
  
})

test_that("bmdHetVar on GiantKelp model, def = hybridAdd, backgtype = hybridSD", {
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  object0 <- drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LL.4())
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  result <- bmdHetVar(object0, bmr = 0.1, backgType = "hybridSD", backg = 2, def = "hybridAdd", R = 50, level = 0.95, progressInfo = FALSE, display = FALSE)
  
  # result
  expect_true(!is.na(result$Results[1, "BMD"]))
  expect_equal(result$Results[1, "BMD"], c(BMD=10.4365491827332))
  expect_equal(result$bmrScaled[1,1], 16.5674974609854)
  expect_equal(result$bmrScaled[1,1], unname(object0$curve(result$Results[1, "BMD"])))
  expect_equal(unname(result$interval[1,]), c(3.4430271147214,21.3097574205103), tolerance = 1e-1)
  
})



