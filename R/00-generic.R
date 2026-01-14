# =============================================================================
# Core Generic Functions
# =============================================================================
#' @import FLCore
#' 
#' Rebuild a fish population
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

