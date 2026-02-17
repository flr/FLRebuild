# Test life history and demographic functions

test_that("covarFn works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate covariates
  covars <- covarFn(brp)
  
  # Should return data.frame or list
  expect_true(is.data.frame(covars) || is.list(covars))
})

test_that("leslieFn works", {
  library(FLCore)
  
  # Create test data
  # This may need specific input format
  # leslie_result <- leslieFn(...)
  # expect_s4_class(leslie_result, "FLPar") or appropriate class
})

test_that("processError works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Add observed data
  ssb.obs(brp) <- ssb(ple4)
  catch.obs(brp) <- catch(ple4)
  
  # Calculate process error
  pe <- processError(brp)
  
  # Should return FLQuant or numeric
  expect_true(is(pe, "FLQuant") || is.numeric(pe) || is.data.frame(pe))
})
