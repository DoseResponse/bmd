# Tests for drmHetVar function
# - Arguments and structure
#   - Missing arguments
#   - Correct backgType and def accepted
# - Ryegrass model
#   - correct values of sigmaFun
#   - correct formula
#   - correct coefficients of sigmaMod
# - GiantKelp model
#   - correct values of sigmaFun
#   - correct formula
#   - correct coefficients of sigmaMod




# Arguments and structure -------------------------------------------------

test_that("drmHetVar handles missing required arguments", {
  expect_error(drmHetVar(), "argument \"formula\" is missing, with no default")
  expect_error(drmHetVar(rootl ~ conc), "argument \"var.formula\" is missing, with no default")
  expect_error(drmHetVar(rootl ~ conc, ~fitted), "argument \"data\" must be supplied")
  expect_error(drmHetVar(rootl ~ conc, ~fitted, data = drcData::ryegrass), "argument \"fct\" must be supplied")
})

test_that("drmHetVar handles missing observations", {
  ryegrass <- drcData::ryegrass
  
  expect_equal(drmHetVar(rootl ~ conc, ~fitted, data = rbind(ryegrass, c(NA, 0)), fct = LL.4())$curvePar,
               c(b = 2.70146939078469, c = 0.411238515917276, d = 7.79533990962011, e = 3.11280193302478), tolerance = 1e-8)
  expect_equal(drmHetVar(rootl ~ conc, ~fitted, data = rbind(ryegrass, c(0, NA)), fct = LL.4())$curvePar,
               c(b = 2.70146939078469, c = 0.411238515917276, d = 7.79533990962011, e = 3.11280193302478), tolerance = 1e-8)
  expect_equal(drmHetVar(rootl ~ conc, ~fitted, data = rbind(ryegrass, c(NA, NA)), fct = LL.4())$curvePar,
               c(b = 2.70146939078469, c = 0.411238515917276, d = 7.79533990962011, e = 3.11280193302478), tolerance = 1e-8)
  expect_equal(drmHetVar(rootl ~ conc, ~fitted, data = cbind(ryegrass, na.col = c(1,NA)), fct = LL.4())$curvePar,
               c(b = 2.70146939078469, c = 0.411238515917276, d = 7.79533990962011, e = 3.11280193302478), tolerance = 1e-8)
  expect_equal(drmHetVar(rootl ~ conc, ~fitted, data = cbind(ryegrass, na.col = NA), fct = LL.4())$curvePar,
               c(b = 2.70146939078469, c = 0.411238515917276, d = 7.79533990962011, e = 3.11280193302478), tolerance = 1e-8)
  
})





# Ryegrass model ----------------------------------------------------------

test_that("drmHetVar on Ryegrass model", {
  var.formula0 <- ~ fitted + I(fitted^2)
  ryegrass.W2.4.hetVar <- drmHetVar(rootl ~ conc, var.formula0, data = drcData::ryegrass, fct = W2.4())
  
  expect_equal(unname(ryegrass.W2.4.hetVar$sigmaFun(0:30)), 
               c(0.445223595938359, 0.452982799562158, 0.684142804788777, 0.729133309006225, 0.62712931196859, 
                 0.517200820094488, 0.429080850573182, 0.362481488223554, 0.312371686102985, 0.274232074339055,
                 0.244734374622287, 0.221538950706607, 0.20300983626506, 0.187992726721281, 0.175661515738903,
                 0.165415590233537, 0.156811341769763, 0.149516137204819, 0.143276993385066, 0.137898950836005,
                 0.133229921433214, 0.129149916870622, 0.125563283907698, 0.12239303259968, 0.119576641433876,
                 0.117062918313419, 0.114809625763844, 0.112781665783847, 0.110949679063262, 0.109288954200143,
                 0.107778571110093), tolerance = 1e-8)
  expect_equal(ryegrass.W2.4.hetVar$var.formula, var.formula0)
  expect_equal(ryegrass.W2.4.hetVar$curvePar,
               c(b = -1.80835016587568, c = 0.251695674718818, d = 7.73421602363691, e = 2.50022750625478), tolerance = 1e-8)
  expect_equal(ryegrass.W2.4.hetVar$sigmaPar, 
               c('(Intercept)' = 0.00730403731249006, 
                 'fitted' = 0.311027994530052, 
                 'I(fitted^2)' = -0.0328936930911057), tolerance = 1e-8)
})



# GiantKelp model ---------------------------------------------------------

test_that("drmHetVar on GiantKelp model", {
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  GiantKelp.LL.4.hetVar <- drmHetVar(tubeLength ~ dose, var.formula0, data = drcData::GiantKelp, fct = LL.4())
  
  expect_equal(unname(GiantKelp.LL.4.hetVar$sigmaFun(0:30)), 
               c(1.11726591406817, 1.55167099990238, 1.73263234952392, 1.82828103973254, 
                 1.88375349592666, 1.91693820464214, 1.93647807446137, 1.94709603355865, 
                 1.95160088283562, 1.95176983121002, 1.94878189680672, 1.94344901375631, 
                 1.93634744430458, 1.92789646121981, 1.91840750107324, 1.90811598138072,
                 1.89720252751968, 1.88580750723786, 1.87404121023516, 1.86199112048934,
                 1.84972720356378, 1.8373058112302, 1.82477260561939, 1.81216477686644,
                 1.79951274423047, 1.78684147457951, 1.77417151400327, 1.7615198019742, 
                 1.74890031900735, 1.73632460564062, 1.72380218110818), tolerance = 1e-8)
  expect_equal(GiantKelp.LL.4.hetVar$var.formula, var.formula0)
  expect_equal(GiantKelp.LL.4.hetVar$curvePar,
               c(b = 1.31919020059578, c = 5.59288746292885, d = 18.2001039179043, e = 44.2426813468012), tolerance = 1e-8)
  expect_equal(GiantKelp.LL.4.hetVar$sigmaPar, 
               c('(Intercept)' = 1.11726591406817, 
                 'log(dose + 1)' = 0.740539128243396, 
                 'I(log(dose + 1)^2)' = -0.164214857054588), tolerance = 1e-8)
})



