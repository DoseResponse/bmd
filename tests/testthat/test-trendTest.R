# Tests for trendTest function
# - Arguments and structure
#   - Missing arguments
# - Simple data 
#   - correct trend test results for all types
# - Ryegrass data
#   - correct trend test results for all types




# Arguments and structure -------------------------------------------------

test_that("trendTest function handles missing required arguments", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(12, 10, 5, 4, 4), 4) + rnorm(20, sd = 0.1))
  
  expect_error(trendTest(x = "x", data = data), 'argument "y" is missing, with no default')
  expect_error(trendTest(y = "y", data = data), 'argument "x" is missing, with no default')
  expect_error(trendTest(x = "x", y = "y", data = data, test = "unknown_test"), "'arg' should be one of \"william\", \"shirley\", \"tukey\"")
})



# Simple data -------------------------------------------------------------

test_that("trendTest (test = 'william') on simple data set with trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(12, 10, 5, 4, 4), 4) + rnorm(20, sd = 0.1))
  result <- trendTest("x", "y", data = data, test = "william")
  
  expect_null(result$p.values)
  expect_equal(result$decisions, matrix(rep("accept", 4), nrow = 4, dimnames = list(c("mu1", "mu2", "mu3", "mu4"), "ctr")))
  expect_true(result$acceptTrend)
})

test_that("trendTest (test = 'william') on simple data set without trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(5, 5, 5, 5, 5), 4) + rnorm(20, sd = 0.1))
  result <- trendTest("x", "y", data = data, test = "william")
  
  expect_null(result$p.values)
  expect_equal(result$decisions, matrix(rep("reject", 4), nrow = 4, dimnames = list(c("mu1", "mu2", "mu3", "mu4"), "ctr")))
  expect_true(!result$acceptTrend)
})

test_that("trendTest (test = 'william') on simple data set with mixed trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(12, 12, 5, 5, 5), 4) + rnorm(20, sd = 0.1))
  result <- trendTest("x", "y", data = data, test = "william")
  
  expect_null(result$p.values)
  expect_equal(result$decisions, matrix(c("reject", "accept", "accept", "accept"), nrow = 4, dimnames = list(c("mu1", "mu2", "mu3", "mu4"), "ctr")))
  expect_true(result$acceptTrend)
})

test_that("trendTest (test = 'shirley') on simple data set with trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(12, 10, 5, 4, 4), 4) + rnorm(20, sd = 0.1))
  result <- trendTest("x", "y", data = data, test = "shirley")
  
  expect_null(result$p.values)
  expect_equal(result$decisions, matrix(c("reject", "accept", "accept", "accept"), nrow = 4, dimnames = list(c("mu1", "mu2", "mu3", "mu4"), "ctr")))
  expect_true(result$acceptTrend)
})

test_that("trendTest (test = 'shirley') on simple data set without trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(5, 5, 5, 5, 5), 4) + rnorm(20, sd = 0.1))
  result <- trendTest("x", "y", data = data, test = "shirley")
  
  expect_null(result$p.values)
  expect_equal(result$decisions, matrix(rep("reject", 4), nrow = 4, dimnames = list(c("mu1", "mu2", "mu3", "mu4"), "ctr")))
  expect_true(!result$acceptTrend)
})

test_that("trendTest (test = 'shirley') on simple data set with mixed trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(12, 12, 5, 5, 5), 4) + rnorm(20, sd = 0.1))
  result <- trendTest("x", "y", data = data, test = "shirley")
  
  expect_null(result$p.values)
  expect_equal(result$decisions, matrix(c("reject", "accept", "accept", "accept"), nrow = 4, dimnames = list(c("mu1", "mu2", "mu3", "mu4"), "ctr")))
  expect_true(result$acceptTrend)
})

test_that("trendTest (test = 'tukey') on simple data set with trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(12, 10, 5, 4, 4), 4) + rnorm(20, sd = 0.1))
  result <- trendTest("x", "y", data = data, test = "tukey")
  
  expect_equal(result$p.values, c("xari: xari" = 0, "xord: xord" = 0, "xarilog: xarilog" = 0))
  expect_equal(result$decisions, c("xari: xari" = "accept", "xord: xord" = "accept", "xarilog: xarilog" = "accept"))
  expect_true(result$acceptTrend)
})

test_that("trendTest (test = 'tukey') on simple data set without trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(5, 5, 5, 5, 5), 4) + rnorm(20, sd = 0.1))
  result <- trendTest("x", "y", data = data, test = "tukey")
  
  expect_equal(result$p.values, c("xari: xari" = 0.634501907632532, "xord: xord" = 0.634521888020178, "xarilog: xarilog" = 0.656543069891262))
  expect_equal(result$decisions, c("xari: xari" = "reject", "xord: xord" = "reject", "xarilog: xarilog" = "reject"))
  expect_true(!result$acceptTrend)
})

test_that("trendTest (test = 'tukey') on simple data set with mixed trend", {
  set.seed(1)
  data <- data.frame(x = rep(1:5, 4), 
                     y = rep(c(12, 12, 5, 5, 5), 4) + rnorm(20, sd = 0.1))
  result <- trendTest("x", "y", data = data, test = "tukey")
  
  expect_equal(result$p.values, c("xari: xari" = 3.67372798848464e-13, "xord: xord" = 3.96238597488718e-13, "xarilog: xarilog" = 1.59872115546023e-14))
  expect_equal(result$decisions, c("xari: xari" = "accept", "xord: xord" = "accept", "xarilog: xarilog" = "accept"))
  expect_true(result$acceptTrend)
})







# Ryegrass ----------------------------------------------------------------

test_that("trendTest (test = 'william') on Ryegrass data set", {
  result <- trendTest("conc", "rootl", data = drcData::ryegrass, test = "william")
  
  expect_null(result$p.values)
  expect_equal(result$decisions, matrix(c("reject", "accept", "accept", "accept", "accept", "accept"), 
                                        nrow = 6, dimnames = list(c("mu1", "mu2", "mu3", "mu4", "mu5", "mu6"), "ctr")))
  expect_true(result$acceptTrend)
})


test_that("trendTest (test = 'shirley') on Ryegrass data set", {
  result <- trendTest("conc", "rootl", data = drcData::ryegrass, test = "shirley")
  
  expect_null(result$p.values)
  expect_equal(result$decisions, matrix(c("reject", "accept", "accept", "accept", "accept", "accept"), 
                                        nrow = 6, dimnames = list(c("mu1", "mu2", "mu3", "mu4", "mu5", "mu6"), "ctr")))
  expect_true(result$acceptTrend)
})

test_that("trendTest (test = 'tukey') on Ryegrass data set", {
  result <- trendTest("conc", "rootl", data = drcData::ryegrass, test = "tukey")
  
  expect_equal(result$p.values, c("xari: xari" = 2.36404196218842e-09, "xord: xord" = 0, "xarilog: xarilog" = 0))
  expect_equal(result$decisions, c("xari: xari" = "accept", "xord: xord" = "accept", "xarilog: xarilog" = "accept"))
  expect_true(result$acceptTrend)
})


