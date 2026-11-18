# S.alba 

test_that("compParmBmd function compares BMD (point) correctly for S.alba model", {
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = S.alba, fct = LL.4())
  
  bmdResult <- bmd(object0, bmr = 3.2, def = "point", backgType = "modelBased", display = FALSE)
  compParmAdditive <- compParmBmd(bmdResult, operator = "-", display = FALSE)
  compParmAdditiveSandwich <- compParmBmd(bmdResult, operator = "-", display = FALSE, vcov. = sandwich::sandwich)
  compParmProportional <- compParmBmd(bmdResult, operator = "/", display = FALSE)
  compParmProportionalSandwich <- compParmBmd(bmdResult, operator = "/", display = FALSE, vcov. = sandwich::sandwich)
  
  # Expected results based on manual calculation (checked in v2.6.8)
  # compParmAdditive
  expect_true(all(!is.na(compParmAdditive)))
  expect_equal(dimnames(compParmAdditive)[[1]], "BMD: Glyphosate-Bentazone")
  expect_equal(dimnames(compParmAdditive)[[2]], c("Estimate","Std. Error","t-value","p-value"))
  expect_equal(as.numeric(compParmAdditive), c(17.314608569836,5.37275063543248,3.2226711687769,0.00205446059829508), tolerance = 1e-6)
  
  # compParmAdditiveSandwich
  expect_true(all(!is.na(compParmAdditiveSandwich)))
  expect_equal(dimnames(compParmAdditiveSandwich)[[1]], "BMD: Glyphosate-Bentazone")
  expect_equal(dimnames(compParmAdditiveSandwich)[[2]], c("Estimate","Std. Error","t-value","p-value"))
  expect_equal(as.numeric(compParmAdditiveSandwich), c(17.314608569836,5.53648925673984,3.12736244340367,0.00272014865341848), tolerance = 1e-6)
  
  # compParmProportional
  expect_true(all(!is.na(compParmProportional)))
  expect_equal(dimnames(compParmProportional)[[1]], "BMD: Glyphosate/Bentazone")
  expect_equal(dimnames(compParmProportional)[[2]], c("Estimate","Std. Error","t-value","p-value"))
  expect_equal(as.numeric(compParmProportional), c(1.7807572610268,0.209102391946808,3.73385141010445,0.000421657961976735), tolerance = 1e-6)
  
  # compParmProportionalSandwich
  expect_true(all(!is.na(compParmProportionalSandwich)))
  expect_equal(dimnames(compParmProportionalSandwich)[[1]], "BMD: Glyphosate/Bentazone")
  expect_equal(dimnames(compParmProportionalSandwich)[[2]], c("Estimate","Std. Error","t-value","p-value"))
  expect_equal(as.numeric(compParmProportionalSandwich), c(1.7807572610268,0.21548358554759,3.62327951357748,0.00060023635604273), tolerance = 1e-6)
  
})

test_that("compParmBmd function compares BMD (relative) correctly for S.alba model", {
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = S.alba, fct = LL.4())
  
  bmdResult <- bmd(object0, bmr = 0.08, def = "relative", backgType = "modelBased", display = FALSE)
  compParmAdditive <- compParmBmd(bmdResult, operator = "-", display = FALSE)
  # Expected results based on manual calculation (checked in v2.6.8)
  # compParmAdditive
  expect_true(all(!is.na(compParmAdditive)))
  expect_equal(dimnames(compParmAdditive)[[1]], "BMD: Glyphosate-Bentazone")
  expect_equal(dimnames(compParmAdditive)[[2]], c("Estimate","Std. Error","t-value","p-value"))
  expect_equal(as.numeric(compParmAdditive), c(9.1055475953203,6.9926792416696,1.3021543360748,0.197840342965524), tolerance = 1e-6)
  
  skip_on_cran()
  set.seed(123)
  bmdBootResult <- bmdBoot(object0, bmr = 0.08, def = "relative", backgType = "modelBased", display = FALSE, R = 100)
  compParmAdditiveBoot <- compParmBmd(bmdBootResult, operator = "-", display = FALSE)
  
  # compParmAdditiveBoot
  expect_true(all(!is.na(compParmAdditiveBoot)[,c(1,2,4)]))
  expect_true(is.na(compParmAdditiveBoot)[,3])
  expect_equal(dimnames(compParmAdditiveBoot)[[1]], "BMD: Glyphosate-Bentazone")
  expect_equal(dimnames(compParmAdditiveBoot)[[2]], c("Estimate","Std. Error","t-value","p-value"))
  expect_equal(as.numeric(compParmAdditiveBoot), c(9.1055475953203,9.16498114107078,NA,0.38), tolerance = 1e-6)
})

test_that("compParmBmd function compares  BMD (hybridExc with hybridSD background) correctly for S.alba model", {
  skip_on_cran()
  object0 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = S.alba, fct = LL.4())
  
  bmdResult <- bmd(object0, bmr = 0.1, def = "hybridExc", backgType = "hybridSD", backg = 2, display = FALSE)
  compParmAdditive <- compParmBmd(bmdResult, operator = "-", display = FALSE)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  # compParmAdditive
  expect_true(all(!is.na(compParmAdditive)))
  expect_equal(dimnames(compParmAdditive)[[1]], "BMD: Glyphosate-Bentazone")
  expect_equal(dimnames(compParmAdditive)[[2]], c("Estimate","Std. Error","t-value","p-value"))
  expect_equal(as.numeric(compParmAdditive), c(8.99619389800992,7.1897490553988,1.25125283632183,0.215699053596006), tolerance = 1e-6)
})


