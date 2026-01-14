# rebuild: Rebuilding Analysis Functions

Consolidated R package containing functions from FLCandy, FLRebuild, and rebuild packages needed for fisheries stock rebuilding analysis workflows.

## Installation

```r
# Install from local source
devtools::install("flr/rebuild")
```

## Purpose

This package consolidates all functions needed for rebuilding analysis, eliminating the need to source multiple files from different packages. All functions previously sourced in analysis scripts are now available through a single package.

## Functions Included

### From FLCandy:
- `toLogits()`, `fromLogits()` - Steepness parameter transformations
- `ftmb()`, `ftmb2()`, `ftmb3()` - Stock-recruitment relationship fitting
- `updateRefs()` - Update reference points
- `pellatParams()`, `pellatParamFn()` - Pella-Thompson production model parameters
- `abi()`, `abiAge()`, `abiMsy()` - Age-based indicators
- `processError()` - Process error calculation
- `covarFn()`, `leslieFn()` - Life history and demographic covariates

### From FLRebuild:
- `invSRR()`, `refCreate()`, `rmax()`, `rmsy()`, `rvirgin()` - Stock-recruitment reference points

### From rebuild:
- Rebuilding trajectory functions (`rebuild()`, `rebuildTime()`, etc.)
- Pella-Tomlinson model functions
- Bootstrap and rebuilding analysis functions

## Usage

```r
library(rebuild)

# All functions are now available without sourcing
library(FLCore)
library(FLBRP)

# Example: Fit stock-recruitment relationship
# sr <- ftmb2(flsr_object, spr0 = 0.7)

# Example: Update reference points
# brp <- updateRefs(brp_object)

# Example: Calculate covariates
# covars <- covarFn(brp_object)
```

## Breaking Changes

⚠️ **Version 0.1.0 includes breaking changes:**
- `from_logits()` → `fromLogits()`
- `to_logits()` → `toLogits()`

See `NEWS.md` for migration guide.

## Package Status

This package consolidates functions from FLCandy and FLRebuild packages. As of version 0.1.0, this package replaces FLRebuild and is self-contained (FLRebuild is no longer a dependency). FLCandy may still be required for some internal methods.

