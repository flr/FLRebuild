# Test calcPriors and getPriors functions

test_that("calcPriors works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Add observed data
  ssb.obs(brp) <- ssb(ple4)
  catch.obs(brp) <- catch(ple4)
  fbar.obs(brp) <- fbar(ple4)
  
  # Calculate priors
  priors <- calcPriors(brp)
  
  # Should return FLPar
  expect_s4_class(priors, "FLPar")
  
  # Should have expected parameters
  expect_true("r" %in% dimnames(priors)$params)
  expect_true("fmsy" %in% dimnames(priors)$params)
  expect_true("bmsy" %in% dimnames(priors)$params)
})

test_that("calcPriors works with FLStock", {
  library(FLCore)
  data(ple4)
  
  # Add benchmark attributes if needed
  # attributes(ple4)$benchmark <- list(Fmsy = 0.2)
  # attributes(ple4)$eqsim <- list(BMSY = 1000, B0 = 2000)
  
  # Calculate priors (may need attributes set)
  # priors <- calcPriors(ple4)
  # expect_s4_class(priors, "FLPar")
})

test_that("getPriors works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Get priors
  priors <- getPriors(brp)
  
  # Should return data.frame
  expect_s3_class(priors, "data.frame")
  
  # Should have expected columns
  expect_true("bmsy" %in% names(priors))
  expect_true("fmsy" %in% names(priors))
  expect_true("msy" %in% names(priors))
})
