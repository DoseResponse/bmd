# Tests for drmHetVar function
# - Arguments and structure
#   - Missing arguments




# Arguments and structure -------------------------------------------------

NULL






# Example -----------------------------------------------------------------
test_that("Example usage of drmMMRE function", {
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  data0 <- data.frame(x = rep(drcData::ryegrass$conc, 2),
                      y = rep(drcData::ryegrass$rootl, 2) +
                        c(rnorm(n = nrow(drcData::ryegrass), mean = 2, sd = 0.5),
                          rnorm(n = nrow(drcData::ryegrass), mean = 2.7, sd = 0.7)),
                      EXP_ID = rep(as.character(1:2), each = nrow(drcData::ryegrass)))
  
  modMMRE <- drmMMRE(y~x, exp_id = EXP_ID, data = data0, fct = LL.4())
  
  expect_equal(coef(modMMRE), 
               c('b:(Intercept)' = 3.14371374389382, 'c:(Intercept)' = 2.96584875840443,
                 'd:(Intercept)' = 10.1745816897618, 'e:(Intercept)' = 2.96457269209798), tolerance = 1e-4)
  expect_equal(vcov(modMMRE), 
               matrix(c(0.419407599823133, -0.0606734434013714, -0.157432089033145, -0.011974063328842,
                        -0.0606734434013714, 0.126805208095318, 0.0802311831922223, -0.0135014412091105, 
                        -0.157432089033145, 0.0802311831922223, 0.134632867305704, -0.00877235238774721, 
                        -0.011974063328842, -0.0135014412091105, -0.00877235238774721, 0.0384415792648907),
                      ncol = 4, nrow = 4,
                      dimnames = list(c("Coefb", "Coefc", "Coefd", "Coefe"),
                                      c("Coefb", "Coefc", "Coefd", "Coefe"))), 
               tolerance = 1e-6)
  
})





