# Test rebuild and rebuildTime functions

test_that("rebuild works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Project rebuilding trajectory
  rebuild_traj <- rebuild(brp,
                          targetF = 0,
                          targetSSB = refpts(brp)["msy", "ssb"])
  
  # Should return FLStock
  expect_s4_class(rebuild_traj, "FLStock")
  
  # Should have SSB
  expect_s4_class(ssb(rebuild_traj), "FLQuant")
  
  # SSB should increase over time (rebuilding)
  ssb_vals <- c(ssb(rebuild_traj))
  if (length(ssb_vals) > 1) {
    # Generally increasing (allowing for some variation)
    expect_true(ssb_vals[length(ssb_vals)] >= ssb_vals[1] * 0.9)
  }
})

test_that("rebuildTime works with FLStock and FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate time to rebuild
  rebuild_time <- rebuildTime(ple4, brp)
  
  # Should return numeric or FLQuant
  expect_true(is.numeric(rebuild_time) || is(rebuild_time, "FLQuant"))
  
  # Should be positive
  if (is.numeric(rebuild_time)) {
    expect_true(rebuild_time > 0 || is.na(rebuild_time))
  }
})

test_that("rebuildTime2 works", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate time to rebuild (method 2)
  rebuild_time <- rebuildTime2(ple4, brp)
  
  # Should return numeric or FLQuant
  expect_true(is.numeric(rebuild_time) || is(rebuild_time, "FLQuant"))
})

test_that("rebuildTime3 works with FLPar", {
  library(FLCore)
  
  # Create FLPar with reference points
  refs <- FLPar(fmsy = 0.2, bmsy = 1000, k = 2000)
  dimnames(refs)$params <- c("fmsy", "bmsy", "k")
  
  # Calculate time to rebuild
  rebuild_time <- rebuildTime3(refs)
  
  # Should return numeric or FLPar
  expect_true(is.numeric(rebuild_time) || is(rebuild_time, "FLPar"))
})

test_that("pBuild works", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate probability of rebuilding
  p_build <- pBuild(ssb(ple4), brp)
  
  # Should return FLQuant or numeric
  expect_true(is(p_build, "FLQuant") || is.numeric(p_build))
  
  # Values should be between 0 and 1 (probabilities)
  if (is(p_build, "FLQuant")) {
    vals <- c(p_build)
    expect_true(all(vals >= 0 & vals <= 1, na.rm = TRUE))
  }
})
