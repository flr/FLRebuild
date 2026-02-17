# Test age-based indicator functions

test_that("abi works with FLStock and FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate ABI
  abi_val <- abi(ple4, brp)
  
  # Should return FLQuant or numeric
  expect_true(is(abi_val, "FLQuant") || is.numeric(abi_val))
  
  # Should have positive values
  if (is(abi_val, "FLQuant")) {
    expect_true(all(abi_val > 0, na.rm = TRUE))
  }
})

test_that("abiAge works with FLStock and FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate ABI by age
  abi_age <- abiAge(ple4, brp)
  
  # Should return FLQuant
  expect_s4_class(abi_age, "FLQuant")
  
  # Should have age dimension
  expect_true(dim(abi_age)[1] > 1)
})

test_that("abiMsy works with FLStock and FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate ABI at MSY
  abi_msy <- abiMsy(ple4, brp)
  
  # Should return FLQuant or numeric
  expect_true(is(abi_msy, "FLQuant") || is.numeric(abi_msy))
})

test_that("forage works with FLStock", {
  library(FLCore)
  data(ple4)
  
  # Calculate forage
  forage_val <- forage(ple4)
  
  # Should return FLQuant
  expect_s4_class(forage_val, "FLQuant")
  
  # Should have positive values
  expect_true(all(forage_val > 0, na.rm = TRUE))
})

test_that("forage works with FLBRP", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Create FLBRP
  brp <- brp(FLBRP(ple4))
  
  # Calculate forage
  forage_val <- forage(brp)
  
  # Should return FLQuant
  expect_s4_class(forage_val, "FLQuant")
})
