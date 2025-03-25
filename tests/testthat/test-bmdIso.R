# Tests for bmdMA function
# - Arguments and structure
#   - Missing arguments
#   - Correct backgType and def accepted
# - formaldehyde model
#   - correct bmd estimate (all definitions for binomial data)
# - ryegrass model
#   - correct bmd estimate (all definitions for cts. data)


# Arguments and structure -------------------------------------------------

## MISSING



# formaldehyde model ------------------------------------------------------

## formaldehyde data example from bmdIso.Rd

test_that("bmdIso function computes BMD (point) correctly for formaldehyde model", {
  formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
                             tumor.incidence = c(0, 0, 0, 3, 21, 150),
                             total = c(122, 27, 126, 113, 34, 182))
  
  result <- bmdIso(tumor.incidence/total ~ conc, 
                   data=formaldehyde, 
                   type="binomial",
                   bmr=0.18,
                   backgType = "modelBased",
                   def = "point")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result))
  expect_equal(result, 7.03841479524439)
})

test_that("bmdIso function computes BMD (excess) correctly for formaldehyde model", {
  formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
                             tumor.incidence = c(0, 0, 0, 3, 21, 150),
                             total = c(122, 27, 126, 113, 34, 182))
  
  result <- bmdIso(tumor.incidence/total ~ conc, 
                   data=formaldehyde, 
                   type="binomial",
                   bmr=0.08,
                   backgType = "modelBased",
                   def = "excess")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result))
  expect_equal(result, 6.36170849845883)
})

test_that("bmdIso function computes BMD (additional) correctly for formaldehyde model", {
  formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
                             tumor.incidence = c(0, 0, 0, 3, 21, 150),
                             total = c(122, 27, 126, 113, 34, 182))
  
  result <- bmdIso(tumor.incidence/total ~ conc, 
                   data=formaldehyde, 
                   type="binomial",
                   bmr=0.08,
                   backgType = "modelBased",
                   def = "additional")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result))
  expect_equal(result, 6.36170849845883)
})




# ryegrass model ----------------------------------------------------------

## ryegrass example from bmdIso.Rd

test_that("bmdIso function computes BMD (point) correctly for ryegrass model", {
  ryegrass1 <- drcData::ryegrass
  ryegrass1$rootl <- 100-ryegrass1$rootl
  
  # Estimating BMD from isotonic regression using relative risk definition and a BMR=0.05
  result <- bmdIso(rootl ~ conc, 
                   data=ryegrass1, 
                   type="continuous",
                   bmr=93,
                   backgType = "modelBased",
                   def = "point")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result))
  expect_equal(result, 1.44279529207128)
})

test_that("bmdIso function computes BMD (relative) correctly for ryegrass model", {
  ryegrass1 <- drcData::ryegrass
  ryegrass1$rootl <- 100-ryegrass1$rootl
  
  # Estimating BMD from isotonic regression using relative risk definition and a BMR=0.05
  result <- bmdIso(rootl ~ conc, 
                   data=ryegrass1, 
                   type="continuous",
                   bmr=0.05,
                   backgType = "modelBased",
                   def = "relative")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result))
  expect_equal(result, 3.68282590840032)
})

test_that("bmdIso function computes BMD (hybridExc with hydridSD background) correctly for ryegrass model", {
  ryegrass1 <- drcData::ryegrass
  ryegrass1$rootl <- 100-ryegrass1$rootl
  
  # Estimating BMD from isotonic regression using relative risk definition and a BMR=0.05
  result <- bmdIso(rootl ~ conc, 
                   data=ryegrass1, 
                   type="continuous",
                   bmr=0.15,
                   backgType = "hybridSD",
                   def = "hybridExc",
                   backg = 2)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(!is.na(result))
  expect_equal(result, 1.14488103599216)
})



