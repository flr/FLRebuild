# Test PellaTomlinson class and methods

test_that("PellaTomlinson object can be created", {
  library(FLCore)
  
  # Create FLPar with parameters
  params <- FLPar(r = 0.3, p = 0.25, virgin = 1000)
  dimnames(params)$params <- c("r", "p", "virgin")
  
  # Create PellaTomlinson object
  pt <- PellaTomlinson(params = params)
  
  # Should be PellaTomlinson object
  expect_s4_class(pt, "PellaTomlinson")
  
  # Should have params slot
  expect_s4_class(params(pt), "FLPar")
})

test_that("msy method works with PellaTomlinson", {
  library(FLCore)
  
  # Create parameters
  params <- FLPar(r = 0.3, p = 0.25, virgin = 1000)
  dimnames(params)$params <- c("r", "p", "virgin")
  
  pt <- PellaTomlinson(params = params)
  
  # Calculate MSY
  msy_val <- msy(pt)
  
  # Should return numeric or FLPar
  expect_true(is.numeric(msy_val) || is(msy_val, "FLPar"))
  
  # Should be positive
  expect_true(all(msy_val > 0, na.rm = TRUE))
})

test_that("bmsy method works with PellaTomlinson", {
  library(FLCore)
  
  # Create parameters
  params <- FLPar(r = 0.3, p = 0.25, virgin = 1000)
  dimnames(params)$params <- c("r", "p", "virgin")
  
  pt <- PellaTomlinson(params = params)
  
  # Calculate BMSY
  bmsy_val <- bmsy(pt)
  
  # Should return numeric or FLPar
  expect_true(is.numeric(bmsy_val) || is(bmsy_val, "FLPar"))
  
  # Should be positive and less than virgin
  expect_true(all(bmsy_val > 0, na.rm = TRUE))
  if (is.numeric(bmsy_val)) {
    expect_true(bmsy_val < params["virgin"])
  }
})

test_that("fmsy method works with PellaTomlinson", {
  library(FLCore)
  
  # Create parameters
  params <- FLPar(r = 0.3, p = 0.25, virgin = 1000)
  dimnames(params)$params <- c("r", "p", "virgin")
  
  pt <- PellaTomlinson(params = params)
  
  # Calculate FMSY
  fmsy_val <- fmsy(pt)
  
  # Should return numeric or FLPar
  expect_true(is.numeric(fmsy_val) || is(fmsy_val, "FLPar"))
  
  # Should be positive
  expect_true(all(fmsy_val > 0, na.rm = TRUE))
})

test_that("production method works with PellaTomlinson", {
  library(FLCore)
  
  # Create parameters
  params <- FLPar(r = 0.3, p = 0.25, virgin = 1000)
  dimnames(params)$params <- c("r", "p", "virgin")
  
  pt <- PellaTomlinson(params = params)
  
  # Create biomass vector
  biomass <- FLQuant(seq(0, 1000, length.out = 100))
  
  # Calculate production
  prod <- production(biomass, pt)
  
  # Should return FLQuant
  expect_s4_class(prod, "FLQuant")
  
  # Should have same dimensions as biomass
  expect_equal(dim(prod), dim(biomass))
})
