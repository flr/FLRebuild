# Test inverse stock-recruitment functions

test_that("invSRR works with FLBRP and FLQuant", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Create recruitment FLQuant
  rec <- FLQuant(seq(100, 1000, length.out = 10))
  
  # Calculate inverse SRR
  ssb_est <- invSRR(brp, rec)
  
  # Should return FLQuant
  expect_s4_class(ssb_est, "FLQuant")
  
  # Should have same dimensions as recruitment
  expect_equal(dim(ssb_est), dim(rec))
})

test_that("invSRR works with FLBRP and FLPar", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Create recruitment FLPar
  rec <- FLPar(rec = 500)
  
  # Calculate inverse SRR
  ssb_est <- invSRR(brp, rec)
  
  # Should return FLPar or FLQuant
  expect_true(is(ssb_est, "FLPar") || is(ssb_est, "FLQuant"))
})

test_that("refCreate works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Create reference point
  new_ref <- refCreate(brp, ref = "custom", value = 1000, quant = "ssb")
  
  # Should return FLPar or modify brp
  expect_true(is(new_ref, "FLPar") || is(new_ref, "FLBRP"))
})
