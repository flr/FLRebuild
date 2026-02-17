# Test time series functions

test_that("tseries works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP with observed data
  brp <- brp(FLBRP(ple4))
  
  # Add observed SSB (required for tseries)
  ssb.obs(brp) <- ssb(ple4)
  catch.obs(brp) <- catch(ple4)
  fbar.obs(brp) <- fbar(ple4)
  
  # Calculate time series
  ts <- tseries(brp)
  
  # Should return FLQuants
  expect_s4_class(ts, "FLQuants")
  
  # Should have expected components
  expect_true("harvest" %in% names(ts))
  expect_true("yield" %in% names(ts))
  expect_true("ssb" %in% names(ts))
})

test_that("tseries works with FLStock", {
  library(FLCore)
  data(ple4)
  
  # Calculate time series
  ts <- tseries(ple4)
  
  # Should return data.frame
  expect_s3_class(ts, "data.frame")
  
  # Should have expected columns
  expect_true("catch" %in% names(ts))
  expect_true("ssb" %in% names(ts))
  expect_true("f" %in% names(ts))
})

test_that("prodPts works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Add observed data
  ssb.obs(brp) <- ssb(ple4)
  catch.obs(brp) <- catch(ple4)
  fbar.obs(brp) <- fbar(ple4)
  
  # Calculate production points
  pp <- prodPts(brp)
  
  # Should return FLQuants
  expect_s4_class(pp, "FLQuants")
})

test_that("prodFn works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate production function
  pf <- prodFn(brp)
  
  # Should return data.frame
  expect_s3_class(pf, "data.frame")
  
  # Should have expected columns
  expect_true("harvest" %in% names(pf))
  expect_true("yield" %in% names(pf))
  expect_true("ssb" %in% names(pf))
})

test_that("MLP works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Add observed data
  ssb.obs(brp) <- ssb(ple4)
  catch.obs(brp) <- catch(ple4)
  fbar.obs(brp) <- fbar(ple4)
  
  # Calculate MLP
  mlp <- MLP(brp)
  
  # Should return data.frame
  expect_s3_class(mlp, "data.frame")
  
  # Should have ssb and yield columns
  expect_true("ssb" %in% names(mlp))
  expect_true("yield" %in% names(mlp))
})
