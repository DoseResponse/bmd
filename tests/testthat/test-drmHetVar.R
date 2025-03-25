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
  lm_object <- lm(y ~ x, 
                  data = data.frame(x = 0:4, 
                                    y = 1:5 + c(-0.4, 0, 0.2, -0.1, 0.13)))
  drm_ryegrass <- drm(rootl ~ conc, data = drcData::ryegrass, fct = LL.4())
  drm_S.alba <- drm(DryMatter ~ Dose, curveid = Herbicide, data = drcData::S.alba, fct = LL.4())
  
  expect_error(drmHetVar(lm_object), 'object must be a dose-response model of class "drc"')
  expect_error(drmHetVar(drm_S.alba), 'dose-response models with multiple curves not supported for heteroscedasticity analysis')
  expect_error(drmHetVar(drm_ryegrass), 'argument "var.formula" is missing, with no default')
  expect_error(drmHetVar(drm_ryegrass, "~ fitted"), 'argument "formula" must be of class "formula"')
})





# Ryegrass model ----------------------------------------------------------

test_that("drmHetVar on Ryegrass model", {
  ryegrass.W2.4 <- drm(rootl ~ conc, data = drcData::ryegrass, fct = W2.4())
  var.formula0 <- ~ fitted + I(fitted^2)
  ryegrass.W2.4.hetVar <- drmHetVar(ryegrass.W2.4, var.formula0)
  
  expect_equal(unname(ryegrass.W2.4.hetVar$sigmaFun(0:30)), 
               c(0.40511221475221, 0.40968626192675, 0.706103004999106, 0.792892779243004, 0.670380066035056, 
                 0.538422954006715, 0.436184506516667, 0.361508041361228, 0.307007516705357, 0.266619836548278,
                 0.236108689714632, 0.212610588112651, 0.194185769216453, 0.179501681344542, 0.167626441952835,
                 0.157896177513736, 0.14982937328875, 0.143070832121776, 0.137354400911807, 0.132477749754008, 
                 0.128285014290574, 0.124654653614345, 0.121490824831162, 0.118717166664911, 0.116272258432784, 
                 0.11410626087336, 0.112178401926649, 0.110455074239762, 0.108908380793708, 0.107515012445684,
                 0.106255373873181), tolerance = 1e-8)
  expect_equal(ryegrass.W2.4.hetVar$var.formula, var.formula0)
  expect_equal(ryegrass.W2.4.hetVar$sigmaMod$coefficients, 
               c('(Intercept)' = -0.0262169518584718, 
                 'fitted' = 0.364206633951366, 
                 'I(fitted^2)' = -0.0399131012240055))
})



# GiantKelp model ---------------------------------------------------------

test_that("drmHetVar on GiantKelp model", {
  GiantKelp.LL.4 <- drm(tubeLength ~ dose, data = drcData::GiantKelp, fct = LL.4())
  var.formula0 <- ~ log(dose+1) + I(log(dose+1)^2)
  GiantKelp.LL.4.hetVar <- drmHetVar(GiantKelp.LL.4, var.formula0)
  
  expect_equal(unname(GiantKelp.LL.4.hetVar$sigmaFun(0:30)), 
               c(1.04778617981885, 1.52105761523325, 1.71774043309804, 1.82140390270584,
                 1.88129829696735, 1.91693199359302, 1.93772576623764, 1.94882504223662, 
                 1.9532883624276, 1.95304997306311, 1.9493924059426, 1.94319840614618,
                 1.93509414356216, 1.92553493643006, 1.91485878842791, 1.90332103382559,
                 1.89111744462622, 1.87840004800314, 1.8652882012118, 1.85187650121705,
                 1.83824053362485, 1.82444111690957, 1.8105274798966, 1.79653967075753, 
                 1.78251040430148, 1.76846649326726, 1.75442996780752, 1.74041895868067,
                 1.72644839956268, 1.71253058960424, 1.69867564707763), tolerance = 1e-8)
  expect_equal(GiantKelp.LL.4.hetVar$var.formula, var.formula0)
  expect_equal(GiantKelp.LL.4.hetVar$sigmaMod$coefficients, 
               c('(Intercept)' = 1.04778617981885, 
                 'log(dose + 1)' = 0.807525479557079, 
                 'I(log(dose + 1)^2)' = -0.179960519480946))
})



