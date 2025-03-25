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

## formaldehyde data example from bmdIsoBoot.Rd

test_that("bmdIsoBoot function computes BMD (point, boot = \"resample\") correctly for formaldehyde model", {
  formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
                             tumor.incidence = c(0, 0, 0, 3, 21, 150),
                             total = c(122, 27, 126, 113, 34, 182))
  set.seed(1)
  result <- bmdIsoBoot(tumor.incidence/total ~ conc, 
                   data=formaldehyde, 
                   type="binomial",
                   bmr=0.18,
                   backgType = "modelBased",
                   def = "point",
                   R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[,"BMD"], 7.03841479524439)
  expect_equal(result[,"BMDL"], 6.80808829160704)
})

test_that("bmdIsoBoot function computes BMD (excess, boot = \"resample\") correctly for formaldehyde model", {
  formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
                             tumor.incidence = c(0, 0, 0, 3, 21, 150),
                             total = c(122, 27, 126, 113, 34, 182))
  set.seed(1)
  result <- bmdIsoBoot(tumor.incidence/total ~ conc, 
                   data=formaldehyde, 
                   type="binomial",
                   bmr=0.08,
                   backgType = "modelBased",
                   def = "excess",
                   R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[, "BMD"], 6.36170849845883)
  expect_equal(result[, "BMDL"], 6.20108949416342)
})

test_that("bmdIsoBoot function computes BMD (additional, boot = \"resample\") correctly for formaldehyde model", {
  formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
                             tumor.incidence = c(0, 0, 0, 3, 21, 150),
                             total = c(122, 27, 126, 113, 34, 182))
  set.seed(1)
  result <- bmdIsoBoot(tumor.incidence/total ~ conc, 
                   data=formaldehyde, 
                   type="binomial",
                   bmr=0.08,
                   backgType = "modelBased",
                   def = "additional",
                   R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[, "BMD"], 6.36170849845883)
  expect_equal(result[, "BMDL"], 6.20108949416342)
})

test_that("bmdIsoBoot function computes BMD (point, boot = \"pseudorandom\") correctly for formaldehyde model", {
  formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
                             tumor.incidence = c(0, 0, 0, 3, 21, 150),
                             total = c(122, 27, 126, 113, 34, 182))
  set.seed(1)
  result <- bmdIsoBoot(tumor.incidence/total ~ conc, 
                       data=formaldehyde, 
                       type="binomial",
                       bmr=0.18,
                       backgType = "modelBased",
                       def = "point",
                       R = 50,
                       boot = "pseudorandom")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[,"BMD"], 7.03841479524439)
  expect_equal(result[,"BMDL"], 6.84192399708765)
})

test_that("bmdIsoBoot function computes BMD (excess, boot = \"pseudorandom\") correctly for formaldehyde model", {
  formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
                             tumor.incidence = c(0, 0, 0, 3, 21, 150),
                             total = c(122, 27, 126, 113, 34, 182))
  set.seed(1)
  result <- bmdIsoBoot(tumor.incidence/total ~ conc, 
                       data=formaldehyde, 
                       type="binomial",
                       bmr=0.08,
                       backgType = "modelBased",
                       def = "excess",
                       R = 50,
                       boot = "pseudorandom")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[, "BMD"], 6.36170849845883)
  expect_equal(result[, "BMDL"], 6.20786436734456)
})

test_that("bmdIsoBoot function computes BMD (additional, boot = \"pseudorandom\") correctly for formaldehyde model", {
  formaldehyde <- data.frame(conc = c(0.0, 0.7, 2.0, 6.0, 10.0, 15.0),
                             tumor.incidence = c(0, 0, 0, 3, 21, 150),
                             total = c(122, 27, 126, 113, 34, 182))
  set.seed(1)
  result <- bmdIsoBoot(tumor.incidence/total ~ conc, 
                       data=formaldehyde, 
                       type="binomial",
                       bmr=0.08,
                       backgType = "modelBased",
                       def = "additional",
                       R = 50,
                       boot = "pseudorandom")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[, "BMD"], 6.36170849845883)
  expect_equal(result[, "BMDL"], 6.20786436734456)
})




# ryegrass model ----------------------------------------------------------

## ryegrass example from bmdIsoBoot.Rd

test_that("bmdIsoBoot function computes BMD (point) correctly for ryegrass model", {
  ryegrass1 <- drcData::ryegrass
  ryegrass1$rootl <- 100-ryegrass1$rootl
  
  # Estimating BMD from isotonic regression using relative risk definition and a BMR=0.05
  set.seed(1)
  result <- bmdIsoBoot(rootl ~ conc, 
                   data=ryegrass1, 
                   type="continuous",
                   bmr=93,
                   backgType = "modelBased",
                   def = "point",
                   R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[, "BMD"], 1.44279529207128)
  expect_equal(result[, "BMDL"], 1.24074240997315)
})

test_that("bmdIsoBoot function computes BMD (relative) correctly for ryegrass model", {
  ryegrass1 <- drcData::ryegrass
  ryegrass1$rootl <- 100-ryegrass1$rootl
  
  # Estimating BMD from isotonic regression using relative risk definition and a BMR=0.05
  set.seed(1)
  result <- bmdIsoBoot(rootl ~ conc, 
                   data=ryegrass1, 
                   type="continuous",
                   bmr=0.05,
                   backgType = "modelBased",
                   def = "relative",
                   R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[, "BMD"], 3.68282590840032)
  expect_equal(result[, "BMDL"], 3.35376062473603)
})

test_that("bmdIsoBoot function computes BMD (hybridExc with hydridSD background) correctly for ryegrass model", {
  ryegrass1 <- drcData::ryegrass
  ryegrass1$rootl <- 100-ryegrass1$rootl
  
  # Estimating BMD from isotonic regression using relative risk definition and a BMR=0.05
  set.seed(1)
  result <- bmdIsoBoot(rootl ~ conc, 
                   data=ryegrass1, 
                   type="continuous",
                   bmr=0.15,
                   backgType = "hybridSD",
                   def = "hybridExc",
                   backg = 2,
                   R = 50)
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[, "BMD"], 1.14488103599216)
  expect_equal(result[, "BMDL"], 0.732577658128982)
})

test_that("bmdIsoBoot function computes BMD (point, boot = \"pseudorandom\") correctly for ryegrass model", {
  ryegrass1 <- drcData::ryegrass
  ryegrass1$rootl <- 100-ryegrass1$rootl
  
  # Estimating BMD from isotonic regression using relative risk definition and a BMR=0.05
  set.seed(1)
  result <- bmdIsoBoot(rootl ~ conc, 
                       data=ryegrass1, 
                       type="continuous",
                       bmr=93,
                       backgType = "modelBased",
                       def = "point",
                       R = 50,
                       boot = "pseudorandom")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[, "BMD"], 1.44279529207128)
  expect_equal(result[, "BMDL"], 0.92019305649066)
})

test_that("bmdIsoBoot function computes BMD (relative, boot = \"pseudorandom\") correctly for ryegrass model", {
  ryegrass1 <- drcData::ryegrass
  ryegrass1$rootl <- 100-ryegrass1$rootl
  
  # Estimating BMD from isotonic regression using relative risk definition and a BMR=0.05
  set.seed(1)
  result <- bmdIsoBoot(rootl ~ conc, 
                       data=ryegrass1, 
                       type="continuous",
                       bmr=0.05,
                       backgType = "modelBased",
                       def = "relative",
                       R = 50,
                       boot = "pseudorandom")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[, "BMD"], 3.68282590840032)
  expect_equal(result[, "BMDL"], 3.12145839108025)
})

test_that("bmdIsoBoot function computes BMD (hybridExc with hydridSD background, boot = \"pseudorandom\") correctly for ryegrass model", {
  ryegrass1 <- drcData::ryegrass
  ryegrass1$rootl <- 100-ryegrass1$rootl
  
  # Estimating BMD from isotonic regression using relative risk definition and a BMR=0.05
  set.seed(1)
  result <- bmdIsoBoot(rootl ~ conc, 
                       data=ryegrass1, 
                       type="continuous",
                       bmr=0.15,
                       backgType = "hybridSD",
                       def = "hybridExc",
                       backg = 2,
                       R = 50,
                       boot = "pseudorandom")
  
  # Expected results based on manual calculation (checked in v2.6.7)
  expect_true(all(!is.na(result)))
  expect_equal(result[, "BMD"], 1.14488103599216)
  expect_equal(result[, "BMDL"], 0.400880981975031)
})



