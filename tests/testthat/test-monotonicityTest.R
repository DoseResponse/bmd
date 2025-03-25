# Tests for monotonicityTest function
# - Arguments and structure
#   - Missing arguments
# - Simple data 
#   - correct monotonicity test results for all types
# - Ryegrass data
#   - correct monotonicity test results for all types




# Arguments and structure -------------------------------------------------

test_that("monotonicityTest function handles missing required arguments", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(12, 10, 5, 4, 4), 4) + rnorm(20, sd = 0.1))
  
  expect_error(monotonicityTest(x = "x", data = data), 'argument "y" is missing, with no default')
  expect_error(monotonicityTest(y = "y", data = data), 'argument "x" is missing, with no default')
  expect_error(monotonicityTest(x = "x", y = "y", data = data, test = "unknown_test"), "'arg' should be one of \"jonckheere\", \"bartholomew\"")
})



# Simple data -------------------------------------------------------------

test_that("monotonicityTest (test = 'bartholomew') on simple data set with trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(12, 10, 5, 4, 4), 4) + rnorm(20, sd = 0.1))
  result <- monotonicityTest("x", "y", data = data, test = "bartholomew")
  
  expect_equal(result$p.value, 0)
  expect_true(result$acceptMonotonicity)
})

test_that("monotonicityTest (test = 'jonckheere') on simple data set with trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(12, 10, 5, 4, 4), 4) + rnorm(20, sd = 0.1))
  result <- monotonicityTest("x", "y", data = data, test = "jonckheere")
  
  expect_equal(result$p.value, 6.21236026663824e-07)
  expect_true(result$acceptMonotonicity)
})

test_that("monotonicityTest (test = 'bartholomew') on simple data set with mixed trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(5, 5, 7, 5, 4), 4) + rnorm(20, sd = 1.1))
  result <- monotonicityTest("x", "y", data = data, test = "bartholomew")
  
  expect_equal(result$p.value, 0.089)
  expect_true(!result$acceptMonotonicity)
})

test_that("monotonicityTest (test = 'jonckheere') on simple data set with mixed trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(5, 5, 7, 5, 4), 4) + rnorm(20, sd = 0.5))
  result <- monotonicityTest("x", "y", data = data, test = "jonckheere")
  
  expect_equal(result$p.value, 0.232501894630627)
  expect_true(!result$acceptMonotonicity)
})


test_that("monotonicityTest (test = 'bartholomew') on simple data set without trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(5, 20) + rnorm(20, sd = 0.5))
  result <- monotonicityTest("x", "y", data = data, test = "bartholomew")
  
  expect_equal(result$p.value, 0.565)
  expect_true(!result$acceptMonotonicity)
})

test_that("monotonicityTest (test = 'jonckheere') on simple data set without trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(5, 20) + rnorm(20, sd = 0.5))
  result <- monotonicityTest("x", "y", data = data, test = "jonckheere")
  
  expect_equal(result$p.value, 0.129415387402834)
  expect_true(!result$acceptMonotonicity)
})





# Ryegrass ----------------------------------------------------------------

test_that("monotonicityTest (test = 'bartholomew') on Ryegrass data set", {
  result <- monotonicityTest("conc", "rootl", data = drcData::ryegrass, test = "bartholomew")
  
  expect_equal(result$p.value, 0)
  expect_true(result$acceptMonotonicity)
})


test_that("monotonicityTest (test = 'jonckheere') on Ryegrass data set", {
  result <- monotonicityTest("conc", "rootl", data = drcData::ryegrass, test = "jonckheere")
  
  expect_equal(result$p.value, 7.16211783824569e-09)
  expect_true(result$acceptMonotonicity)
})

