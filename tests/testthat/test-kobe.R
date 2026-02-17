# Test Kobe plot functions

test_that("kobeSummary works with data.frame", {
  library(FLCore)
  
  # Create test data
  stock_data <- data.frame(
    year = 2000:2010,
    ssb = seq(0.5, 1.5, length.out = 11),
    harvest = seq(0.3, 1.2, length.out = 11)
  )
  
  # Calculate Kobe summary
  kobe_sum <- kobeSummary(stock_data)
  
  # Should return data.frame
  expect_s3_class(kobe_sum, "data.frame")
  
  # Should have expected columns
  expect_true("ssb" %in% names(kobe_sum) || "stock" %in% names(kobe_sum))
  expect_true("harvest" %in% names(kobe_sum) || "f" %in% names(kobe_sum))
})

test_that("prob works with numeric vectors", {
  library(FLCore)
  
  # Create test data
  stock <- c(0.5, 0.7, 0.9, 1.1, 1.3)
  harvest <- c(0.3, 0.5, 0.7, 0.9, 1.1)
  
  # Calculate probabilities
  prob_vals <- prob(stock, harvest)
  
  # Should return numeric
  expect_type(prob_vals, "double")
  
  # Should have same length as input
  expect_length(prob_vals, length(stock))
  
  # Values should be between 0 and 1 (probabilities)
  expect_true(all(prob_vals >= 0 & prob_vals <= 1, na.rm = TRUE))
})

test_that("results works with FLStock and FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate results
  res <- results(ple4, brp)
  
  # Should return FLQuants
  expect_s4_class(res, "FLQuants")
  
  # Should have expected components
  expect_true("SSB" %in% names(res))
  expect_true("F" %in% names(res))
  expect_true("Yield" %in% names(res))
})
