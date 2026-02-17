# Test natural mortality functions

test_that("M1Fn works with FLQuant", {
  library(FLCore)
  
  # Create FLQuant with natural mortality
  m <- FLQuant(0.2, dimnames = list(age = 1:10))
  
  # Calculate M1
  m1 <- M1Fn(m)
  
  # Should return FLQuant
  expect_s4_class(m1, "FLQuant")
  
  # Should have same dimensions
  expect_equal(dim(m1), dim(m))
  
  # Should be positive
  expect_true(all(m1 > 0, na.rm = TRUE))
})

test_that("M1Fn works with FLStock", {
  library(FLCore)
  data(ple4)
  
  # Calculate M1
  m1 <- M1Fn(ple4)
  
  # Should return FLQuant
  expect_s4_class(m1, "FLQuant")
})

test_that("M2Fn works with FLQuant", {
  library(FLCore)
  
  # Create FLQuant with natural mortality
  m <- FLQuant(0.2, dimnames = list(age = 1:10))
  
  # Calculate M2
  m2 <- M2Fn(m)
  
  # Should return FLQuant
  expect_s4_class(m2, "FLQuant")
  
  # Should have same dimensions
  expect_equal(dim(m2), dim(m))
})

test_that("ddM works with FLQuant and FLPar", {
  library(FLCore)
  
  # Create weight FLQuant
  wt <- FLQuant(seq(0.1, 2, length.out = 10), dimnames = list(age = 1:10))
  
  # Create parameters
  par <- FLPar(m1 = 0.3, m2 = -0.3)
  dimnames(par)$params <- c("m1", "m2")
  
  # Calculate density-dependent mortality
  dd_m <- ddM(wt, par)
  
  # Should return FLQuant
  expect_s4_class(dd_m, "FLQuant")
  
  # Should have same dimensions as weight
  expect_equal(dim(dd_m), dim(wt))
  
  # Should be positive
  expect_true(all(dd_m > 0, na.rm = TRUE))
})
