# =============================================================================
# Core Generic Functions
# =============================================================================
#' Import FLCore but keep package-local S4 generics for names whose FLCore
#' signatures differ (r, production, msy, fmsy, bmsy, fwd, p).
#' @rawNamespace import(FLCore, except = c(r, production, msy, fmsy, bmsy, fwd, p))
#' 
#' @title Rebuild Trajectories
#' 
#' @description Projects rebuilding trajectories from different initial SSB levels
#'
#' @param object An object representing the population
#' @param ... Additional arguments
#' @return An object with rebuilding trajectories
#' @export
setGeneric("rebuild", function(object, ...) standardGeneric("rebuild"))

#' Calculate rebuilding time
#'
#' @param object An object containing rebuilding trajectories (FLStock, biodyn, or FLPar)
#' @param ... Additional arguments
#' @return A data frame with columns year and initial (or T_recover and initial)
#' 
#' @examples
#' \dontrun{
#' # Example with FLStock object
#' library(FLCore)
#' data(ple4)
#' rebuild_times <- rebuildTime(ple4)
#' 
#' # Example with biodyn object
#' library(mpb)
#' params <- FLPar(r=0.5, k=1000, p=1)
#' bd <- biodyn(params)
#' rebuild_times <- rebuildTime(bd)
#' }
#' 
#' @export
setGeneric("rebuildTime", function(object, ...) standardGeneric("rebuildTime"))

#' Calculate Rebuild Time using FLPar parameters (version 2)
#'
#' @param object An object with parameters
#' @param ... Additional arguments
#' @return A data frame with rebuild times
#' @export
setGeneric("rebuildTime2", function(object, ...) standardGeneric("rebuildTime2"))

#' Calculate Rebuild Time using FLPar parameters (version 3)
#'
#' @param object An object with parameters
#' @param ... Additional arguments
#' @return A data frame with rebuild times
#' @export
setGeneric("rebuildTime3", function(object, ...) standardGeneric("rebuildTime3"))

# =============================================================================
# Age-Based Indicator Generic Functions
# =============================================================================

#' Calculate the reference age for a FLBRP object
#'
#' @description This function calculates the reference age for a FLBRP object.
#'
#' @param object FLBRP object.
#' @param ref Reference point, e.g., "msy" (default) or "f0.1".
#' @param p Probability threshold (default = 0.9).
#' @param ... Additional arguments
#'
#' @return An FLQuant object containing reference ages.
#' @export
setGeneric("abiAge", function(object, ref = "msy", p = 0.9, ...) {
  standardGeneric("abiAge")
})

#' Calculate P(N) at FMSY for a FLBRP object
#'
#' @description This function calculates P(N) at FMSY for a FLBRP object.
#'
#' @param object A FLBRP object.
#' @param ref Reference point, e.g., "msy" (default) or "f0.1".
#' @param p Probability threshold (default = 0.9).
#' @param ... Additional arguments
#'
#' @return An FLQuant object containing P(N) at FMSY.
#' @export
setGeneric("abiMsy", function(object, ref = "msy", p = 0.9, ...) {
  standardGeneric("abiMsy")
})

#' Calculate P obs for a FLStock object
#'
#' @description This function calculates P obs for a FLStock object.
#'
#' @param object A FLStock object.
#' @param age Reference ages obtained from abiAge (FLBRP or FLQuant object).
#' @param ... Additional arguments
#'
#' @return An FLQuant object containing P obs.
#' @export
setGeneric("abi", function(object, age, ...) {
  standardGeneric("abi")
})

# =============================================================================
# Process Error Generic Function
# =============================================================================

#' Calculate process error
#'
#' @description This function calculates process error for an FLBRP object.
#'
#' @param object An FLBRP object.
#' @param ... Additional arguments
#'
#' @return An FLQuants object containing process error metrics.
#' @export
setGeneric("processError", function(object, ...) {
  standardGeneric("processError")
})

#' Process error from stock and surplus production
#'
#' Shared S4 generic. The \code{FLStock} methods live in this package
#' (and \code{icesdata}); SS3 methods live in \code{sstruth}.
#'
#' For Pella–Tomlinson dynamics the process error can be placed after catch,
#' before catch, or on productivity (see \code{\link{pe,FLStock,FLBRP-method}}).
#' Catch \eqn{C_y} is used in place of \eqn{F_y B_y}.
#'
#' @param object An \code{FLStock}, SS3 directory, or \code{curveSS} list
#' @param eq An \code{FLBRP}, \code{FLPar} (PT params), or \code{PellaTomlinson};
#'   omitted for SS3 methods
#' @param ... Additional arguments (\code{location}, \code{scale},
#'   \code{production} (\code{"FLBRP"} or \code{"pellat"}), \code{stock})
#' @export
setGeneric("pe", function(object, eq, ...) standardGeneric("pe"))

# =============================================================================
# Forward Projection Generic Function
# =============================================================================

#' Forward project a population model
#'
#' @description Forward projects a population model (e.g., PellaTomlinson) over time.
#'
#' @param object A population model object (e.g., PellaTomlinson).
#' @param ... Additional arguments (F_val, catch, years, initial_biomass, etc.)
#'
#' @return Projected population trajectory.
#' @export
setGeneric("fwd", function(object, ...) {
  standardGeneric("fwd")
})

# =============================================================================
# Analysis Results Generic Function
# =============================================================================

#' Calculate Performance Metrics from Operating Model Results
#'
#' @description
#' Calculates a suite of performance metrics from operating model (OM) results
#' relative to reference points from an equilibrium (eq) object.
#'
#' @param object An FLStock object or list of FLStock objects (operating model results)
#' @param eq An FLBRP object containing reference points
#' @param ... Additional arguments
#'
#' @return An FLQuants object containing relative metrics
#'
#' @export
setGeneric("results", function(object, eq, ...) standardGeneric("results"))

# =============================================================================
# FLBRP Reference Point Generic Functions
# =============================================================================

#' Calculate MSY and Virgin State Metrics
#'
#' @description 
#' Calculates key metrics (F, SSB, catch, ebiomass) at virgin and MSY states for an FLBRP object
#'
#' @param object An object of class FLBRP
#'
#' @return A named vector containing metrics for virgin and MSY states:
#'   \item{virgin.f}{F at virgin state}
#'   \item{virgin.catch}{Catch at virgin state}
#'   \item{virgin.ebiomass}{Equilibrium biomass at virgin state}
#'   \item{msy.f}{F at MSY}
#'   \item{msy.catch}{Catch at MSY}
#'   \item{msy.ebiomass}{Equilibrium biomass at MSY}
#'
#' @note The function excludes virgin.ssb and msy.ssb from the output
#'
#' @export
#'
#' @rdname msyVirgin
#' 
#' @examples
#' \dontrun{
#' data(ple4brp)
#' msyVirgin(ple4brp)
#' }
setGeneric("msyVirgin", function(object, ...) standardGeneric("msyVirgin"))

#' Calculate Blim Reference Point
#' 
#' @param object An FLBRP object
#' @param ratio Ratio of virgin recruitment to use (default 0.3)
#' @return An FLPar object containing Blim reference points
#' @export
setGeneric("blim", function(object, ...) standardGeneric("blim"))

#' Extract Time Series Statistics from FLStock Objects
#'
#' @description
#' Creates a data frame of time series statistics from FLStock objects, including
#' catch, ebiomass, SSB, fishing mortality, harvest rate, and mean natural mortality.
#'
#' @param object An object of class FLStock or FLStocks
#' @param ... Additional arguments (not currently used)
#'
#' @return A data frame containing time series of:
#'   \itemize{
#'     \item catch - Catch values
#'     \item eb - Exploitable biomass
#'     \item ssb - Spawning stock biomass
#'     \item f - Fishing mortality (Fbar)
#'     \item h - Harvest rate (catch/ebiomass)
#'     \item m - Mean natural mortality
#'   }
#'
#' @examples
#' \dontrun{
#' # For single stock
#' data(ple4)
#' ts1 <- tseries(ple4)
#'
#' # For multiple stocks
#' stocks <- FLStocks(stock1=ple4, stock2=ple4)
#' ts2 <- tseries(stocks)
#' }
#'
#' @export
setGeneric("tseries", function(object, ...) standardGeneric("tseries"))

# =============================================================================
# Stock Synthesis Helper Generics
# =============================================================================

#' Run Stock Synthesis diagnostics
#'
#' @param object A Stock Synthesis object or path-like input, depending on method.
#' @param ... Additional arguments passed to methods.
#' @return Method-specific diagnostics output.
#' @export
setGeneric("ssDiagnostics", function(object, ...) standardGeneric("ssDiagnostics"))

#' Run Stock Synthesis retrospective diagnostics
#'
#' @param object A retrospective object (typically a list of SS outputs).
#' @param ... Additional arguments passed to methods.
#' @return Method-specific retrospective diagnostics output.
#' @export
setGeneric("ssRetrospective", function(object, ...) standardGeneric("ssRetrospective"))

#' Build equilibrium and production-curve inputs from SS output
#'
#' @param object A Stock Synthesis output object.
#' @param ... Additional arguments passed to methods.
#' @return Method-specific equilibrium/production-curve output.
#' @export
setGeneric("curveSS", function(object, ...) standardGeneric("curveSS"))

#' Estimate vulnerable biomass from SS report components
#'
#' @param object A Stock Synthesis output object.
#' @param ... Additional arguments passed to methods.
#' @return Method-specific vulnerable biomass output.
#' @export
setGeneric("vBio", function(object, ...) standardGeneric("vBio"))

#' Standardize SS index series and fit smooth trend
#'
#' @param object A data.frame containing index observations.
#' @param ... Additional arguments passed to methods.
#' @return Method-specific standardized index output.
#' @export
setGeneric("ssIndexScale", function(object, ...) standardGeneric("ssIndexScale"))

#' Compute runs-test diagnostics for standardized SS index residuals
#'
#' @param object A data.frame containing standardized index residuals.
#' @param ... Additional arguments passed to methods.
#' @return Method-specific runs-test output.
#' @export
setGeneric("ssIndexRuns", function(object, ...) standardGeneric("ssIndexRuns"))

#' Build plotting data for standardized SS indices
#'
#' @param object A data.frame containing index observations.
#' @param ... Additional arguments passed to methods.
#' @return Method-specific plotting-ready data.
#' @export
setGeneric("ssIndexDiagnostics", function(object, ...) standardGeneric("ssIndexDiagnostics"))

#' Prepare SS index data for diagnostics
#'
#' @param object Raw SS index object (e.g., list output from cpueSS) or a data.frame.
#' @param ... Additional arguments passed to methods.
#' @return Method-specific standardized index input data.
#' @export
setGeneric("ssIndexPrep", function(object, ...) standardGeneric("ssIndexPrep"))

# =============================================================================
# ICES stock-attribute accessors (eqsim, benchmark, fishlife)
# Also exported by icesdata. Defined here so FLRebuild does not need FLCandy.
# =============================================================================

#' Extract EqSim reference points from an FLStock
#'
#' @param object An \code{FLStock}, \code{FLStocks}, or list of stocks
#' @param ... Unused
#' @return An \code{FLPar} (single stock) or a data.frame (several stocks)
#' @export
setGeneric("eqsim", function(object, ...) standardGeneric("eqsim"))

#' Extract ICES benchmark reference points from an FLStock
#'
#' @param object An \code{FLStock}, \code{FLStocks}, \code{FLBRP}, or list
#' @param ... Unused
#' @return An \code{FLPar} (single stock) or a data.frame (several stocks)
#' @export
setGeneric("benchmark", function(object, ...) standardGeneric("benchmark"))

#' Extract FishLife parameters from an FLStock
#'
#' @param object An \code{FLStock}, \code{FLStocks}, or list of stocks
#' @param ... Unused
#' @return An \code{FLPar} (single stock) or a data.frame (several stocks)
#' @export
setGeneric("fishlife", function(object, ...) standardGeneric("fishlife"))

# =============================================================================
# Assessment / projection run loaders
# =============================================================================

#' Load assessment or projection runs
#'
#' S4 generic for named lists of Stock Synthesis projection outputs or JABBA
#' assessment fits. \code{character} methods take a directory or a named vector
#' of paths; \code{list} methods validate an already-loaded collection.
#' Source is chosen with \code{source} (\code{"ss"}, \code{"jabba"}, or
#' \code{"auto"}). Auto-detection uses SS catch-level names / \code{*_\\{catch\\}t.Rdata}
#' files versus JABBA run labels / \code{*_jabba.rdata} files.
#'
#' @param object Named character vector of paths, a directory, or a named list
#'   of loaded run objects.
#' @param ... Passed to methods (\code{pattern}, \code{runs}, \code{source}).
#' @return A named list of run objects (SS_output-shaped or JABBA fits).
#' @examples
#' \dontrun{
#' runs <- getRuns(c(
#'   "0" = "path/to/HW2e_..._0t.Rdata",
#'   "250" = "path/to/HW2e_..._250t.Rdata"))
#' runs <- getRuns("P:/.../SSoutput")
#' jb   <- getRuns("P:/.../outputCorrected", source = "jabba",
#'                 runs = c("1-B", "1-S", "2-S"))
#' }
#' @seealso \code{\link{getJabbaRuns}}, \code{\link{kobeMomentGrid}},
#'   \code{\link{kobeMomentGridJabba}}, \code{\link{makeK2SM}}
#' @export
setGeneric("getRuns", function(object, ...) standardGeneric("getRuns"))
