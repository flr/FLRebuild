# FLRebuild Unit Tests

This directory contains unit tests for the FLRebuild package using the `testthat` framework.

## Running Tests

To run all tests:

```r
library(testthat)
library(FLRebuild)
test_dir("tests/testthat")
```

Or using devtools:

```r
devtools::test()
```

## Test Files

### Core Functionality
- `test-helpers.R` - Helper functions (fromLogits, toLogits, interp, tryIt)
- `test-ebiomass.R` - Exploitable biomass calculations
- `test-spr0Yr.R` - Spawning biomass per recruit calculations
- `test-eql.R` - Equilibrium model fitting

### Reference Points
- `test-refpts.R` - Reference point calculations (rmax, rmsy, rvirgin, msyVirgin, blim, refptsEB)
- `test-calcPriors.R` - Prior calculations for production models

### Rebuilding Analysis
- `test-rebuild.R` - Rebuilding trajectory projection and time calculations
- `test-tseries.R` - Time series and production analysis

### Production Models
- `test-pellatParams.R` - Pella-Tomlinson parameter estimation
- `test-pellaTomlinson.R` - PellaTomlinson class and methods

### Age-Based Indicators
- `test-abi.R` - Age-based indicator calculations

### Supporting Functions
- `test-mortality.R` - Natural mortality functions (M1Fn, M2Fn, ddM)
- `test-kobe.R` - Kobe plot and status assessment
- `test-covariates.R` - Life history covariates
- `test-invSRR.R` - Inverse stock-recruitment functions

## Test Coverage

The tests cover:
- Input validation
- Output type checking
- Dimension consistency
- Value range checks (e.g., positive values, probabilities between 0-1)
- Basic functionality of main methods

## Adding New Tests

When adding new functionality:
1. Create a new test file following the naming convention `test-<functionality>.R`
2. Use descriptive test names with `test_that()`
3. Test both success cases and edge cases
4. Include input validation tests
5. Test with example data (e.g., `ple4` from FLCore)

## Notes

- Some tests may require specific data attributes or setup
- Tests using TMB compilation may take longer
- Some functions may need mock data or specific FLR object configurations
- Tests are designed to be run with example data from FLCore package
