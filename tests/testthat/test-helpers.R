# Test helper functions

test_that("fromLogits converts logit values correctly", {
  # Test basic conversion
  logit_vals <- c(0, 1, -1, 2, -2)
  result <- fromLogits(logit_vals)
  
  # Should return values between 0.2001 and 0.9999
  expect_true(all(result >= 0.2001 & result <= 0.9999))
  
  # Test edge cases
  expect_equal(fromLogits(0), 0.2001 + 0.7998 * 0.5, tolerance = 1e-6)
  
  # Test that large positive values approach 0.9999
  expect_true(fromLogits(10) > 0.99)
  
  # Test that large negative values approach 0.2001
  expect_true(fromLogits(-10) < 0.3)
})

test_that("toLogits converts steepness values correctly", {
  # Test basic conversion
  steepness_vals <- c(0.5, 0.7, 0.9)
  result <- toLogits(steepness_vals)
  
  # Should return numeric values
  expect_type(result, "double")
  expect_length(result, 3)
  
  # Test that conversion is approximately invertible
  original <- c(0.5, 0.7, 0.9)
  logit_vals <- toLogits(original)
  back_converted <- fromLogits(logit_vals)
  expect_equal(back_converted, original, tolerance = 1e-4)
  
  # Test edge cases
  expect_true(is.finite(toLogits(0.2)))
  expect_true(is.finite(toLogits(1.0)))
})

test_that("fromLogits and toLogits are inverse functions", {
  # Test round-trip conversion
  original <- seq(0.21, 0.99, by = 0.1)
  logit_vals <- toLogits(original)
  converted_back <- fromLogits(logit_vals)
  
  expect_equal(converted_back, original, tolerance = 1e-4)
})

test_that("interp function works correctly", {
  # Create test data
  df <- data.frame(
    initial = rep(c(0.1, 0.2, 0.3), each = 5),
    year = rep(1:5, 3),
    biomass = c(0.5, 0.8, 1.2, 1.5, 1.8,  # initial = 0.1, reaches 1 at year 3
                0.3, 0.6, 0.9, 1.1, 1.4,  # initial = 0.2, reaches 1 at year 4
                0.2, 0.4, 0.7, 0.95, 1.2) # initial = 0.3, reaches 1 at year 5
  )
  
  result <- interp(df)
  
  # Should return data frame with initial and year columns
  expect_s3_class(result, "data.frame")
  expect_true(all(c("initial", "year") %in% names(result)))
  
  # Check that it finds the correct years
  expect_equal(result$year[result$initial == 0.1], 3)
  expect_equal(result$year[result$initial == 0.2], 4)
  expect_equal(result$year[result$initial == 0.3], 5)
})

test_that("interp handles cases where biomass never reaches 1", {
  # Create test data where biomass never reaches 1
  df <- data.frame(
    initial = rep(c(0.1, 0.2), each = 3),
    year = rep(1:3, 2),
    biomass = c(0.5, 0.6, 0.7,  # initial = 0.1, never reaches 1
                0.3, 0.4, 0.5)  # initial = 0.2, never reaches 1
  )
  
  result <- interp(df)
  
  # Should return NA for years where biomass never reaches 1
  expect_true(is.na(result$year[result$initial == 0.1]))
  expect_true(is.na(result$year[result$initial == 0.2]))
})

test_that("tryIt function handles errors gracefully", {
  # Test with valid input
  expect_equal(tryIt(5), 5)
  expect_equal(tryIt(c(1, 2, 3)), c(1, 2, 3))
  
  # Test with expression that might fail
  result <- tryIt(log(-1))
  # Should return NA or handle error gracefully
  expect_true(is.na(result) || is.null(result) || inherits(result, "try-error"))
})
