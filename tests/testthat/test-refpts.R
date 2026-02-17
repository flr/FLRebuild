# Test reference point functions

test_that("refptsEB works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate reference points with exploitable biomass
  refpts_eb <- refptsEB(brp)
  
  # Should return FLPar
  expect_s4_class(refpts_eb, "FLPar")
  
  # Should have harvest and yield columns
  expect_true("harvest" %in% dimnames(refpts_eb)$quant)
  expect_true("yield" %in% dimnames(refpts_eb)$quant)
})

test_that("rmax works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate rmax
  rmax_val <- rmax(brp)
  
  # Should return numeric or FLPar
  expect_true(is.numeric(rmax_val) || is(rmax_val, "FLPar"))
  
  # Should be positive
  expect_true(all(rmax_val > 0, na.rm = TRUE))
})

test_that("rmsy works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate rmsy
  rmsy_val <- rmsy(brp)
  
  # Should return numeric or FLPar
  expect_true(is.numeric(rmsy_val) || is(rmsy_val, "FLPar"))
  
  # Should be positive
  expect_true(all(rmsy_val > 0, na.rm = TRUE))
})

test_that("rvirgin works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate rvirgin
  rvirgin_val <- rvirgin(brp)
  
  # Should return numeric or FLPar
  expect_true(is.numeric(rvirgin_val) || is(rvirgin_val, "FLPar"))
  
  # Should be positive
  expect_true(all(rvirgin_val > 0, na.rm = TRUE))
})

test_that("msyVirgin works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate msyVirgin
  msy_virgin <- msyVirgin(brp)
  
  # Should return FLPar or similar
  expect_true(is(msy_virgin, "FLPar") || is.numeric(msy_virgin))
})

test_that("blim works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate blim
  blim_val <- blim(brp)
  
  # Should return numeric or FLPar
  expect_true(is.numeric(blim_val) || is(blim_val, "FLPar"))
  
  # Should be positive
  expect_true(all(blim_val > 0, na.rm = TRUE))
})
