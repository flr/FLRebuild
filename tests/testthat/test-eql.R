# Test eql (equilibrium) function

test_that("eql works with FLStock and bevholtSV model", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Fit equilibrium model
  eql_result <- eql(ple4, model = "bevholtSV")
  
  # Should return FLBRP object
  expect_s4_class(eql_result, "FLBRP")
  
  # Should have reference points
  expect_s4_class(refpts(eql_result), "FLPar")
  
  # Should have MSY reference point
  expect_true("msy" %in% dimnames(refpts(eql_result))$refpt)
  
  # Should have attributes
  expect_true("sr" %in% names(attributes(eql_result)))
  expect_true("logLik" %in% names(attributes(eql_result)))
})

test_that("eql works with rickerSV model", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Fit equilibrium model with Ricker
  eql_result <- eql(ple4, model = "rickerSV")
  
  # Should return FLBRP object
  expect_s4_class(eql_result, "FLBRP")
  
  # Should have reference points
  expect_s4_class(refpts(eql_result), "FLPar")
})

test_that("eql works with segreg model", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Fit equilibrium model with Segreg
  eql_result <- eql(ple4, model = "segreg")
  
  # Should return FLBRP object
  expect_s4_class(eql_result, "FLBRP")
  
  # Should have reference points
  expect_s4_class(refpts(eql_result), "FLPar")
})

test_that("eql accepts prior parameters", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  
  # Fit with priors
  eql_result <- eql(ple4, 
                    model = "bevholtSV",
                    prior_s = 0.7,
                    cv_s = 0.1,
                    prior_r0 = 1000,
                    cv_r0 = 0.2)
  
  # Should return FLBRP object
  expect_s4_class(eql_result, "FLBRP")
})
