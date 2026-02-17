# Test spr0Yr function

test_that("spr0Yr works with FLStock", {
  library(FLCore)
  data(ple4)
  
  # Calculate SPR0
  spr0 <- spr0Yr(ple4)
  
  # Should return FLQuant
  expect_s4_class(spr0, "FLQuant")
  
  # Should have positive values
  expect_true(all(spr0 > 0, na.rm = TRUE))
  
  # Should have same year dimension as stock
  expect_equal(dim(spr0)[2], dim(ple4)[2])
})
