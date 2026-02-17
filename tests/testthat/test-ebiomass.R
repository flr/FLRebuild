# Test exploitable biomass functions

test_that("ebiomass works with FLStock", {
  library(FLCore)
  data(ple4)
  
  # Calculate exploitable biomass
  eb <- ebiomass(ple4)
  
  # Should return FLQuant
  expect_s4_class(eb, "FLQuant")
  
  # Should have correct dimensions (no age dimension)
  expect_equal(length(dim(eb)), 6)
  expect_equal(dim(eb)[1], 1)  # No age dimension
  
  # Should have positive values
  expect_true(all(eb > 0, na.rm = TRUE))
})

test_that("ebiomass works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP object
  brp <- brp(FLBRP(ple4))
  
  # Calculate exploitable biomass
  eb <- ebiomass(brp)
  
  # Should return FLQuant
  expect_s4_class(eb, "FLQuant")
  
  # Should have positive values
  expect_true(all(eb > 0, na.rm = TRUE))
})
