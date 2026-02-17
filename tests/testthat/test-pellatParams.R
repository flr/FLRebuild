# Test Pella-Tomlinson parameter functions

test_that("pellatParams works with FLPar", {
  library(FLCore)
  
  # Create FLPar with reference points
  refs <- FLPar(fmsy = 0.2, bmsy = 1000, k = 2000)
  dimnames(refs)$params <- c("fmsy", "bmsy", "k")
  
  # Calculate parameters
  params <- pellatParams(refs)
  
  # Should return FLPar
  expect_s4_class(params, "FLPar")
  
  # Should have r, k, p parameters
  expect_true("r" %in% dimnames(params)$params)
  expect_true("k" %in% dimnames(params)$params)
  expect_true("p" %in% dimnames(params)$params)
  
  # Parameters should be positive
  expect_true(all(params > 0, na.rm = TRUE))
})

test_that("pellatParams works with numeric vector", {
  # Create named numeric vector
  refs <- c(fmsy = 0.2, bmsy = 1000, virgin = 2000)
  
  # Calculate parameters
  params <- pellatParams(refs)
  
  # Should return FLPar
  expect_s4_class(params, "FLPar")
  
  # Should have required parameters
  expect_true("r" %in% dimnames(params)$params)
  expect_true("k" %in% dimnames(params)$params)
  expect_true("p" %in% dimnames(params)$params)
})

test_that("pellatParams works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP and calculate reference points
  brp <- brp(FLBRP(ple4))
  
  # Calculate parameters
  params <- pellatParams(brp)
  
  # Should return FLPar
  expect_s4_class(params, "FLPar")
  
  # Should have required parameters
  expect_true("r" %in% dimnames(params)$params || "k" %in% dimnames(params)$params)
})
