# =============================================================================
# Mortality Functions (M Functions)
# =============================================================================

#' Calculate Density-Dependent Mortality
#'
#' @description
#' Calculates density-dependent mortality rate based on weight and parameters.
#' Uses the formula: M = m1 * (wt^m2), where m1 and m2 are parameters.
#'
#' @param wt An FLQuant object containing weights
#' @param par An FLPar object containing parameters 'm1' and 'm2'
#' @return An FLQuant object with density-dependent mortality rates
#'
#' @details
#' This function calculates mortality as a power function of weight:
#' M = m1 * (wt^m2)
#' where m1 is the scaling parameter and m2 is the exponent.
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' wt=FLQuant(c(0.1, 0.5, 1.0, 2.0))
#' par=FLPar(m1 = 0.3, m2 = -0.3)
#' ddM(wt, par)
#' }
#'
#' @export
setGeneric("ddM", function(wt, par) standardGeneric("ddM"))

#' @rdname ddM
#' @export
setMethod("ddM", signature(wt = "FLQuant", par = "FLPar"),
  function(wt, par) {
    par["m1"] %*% (wt %^% par["m2"])
  })

#' Calculate Background Mortality (M1)
#'
#' @description
#' Calculates or extracts the background/constant component of natural mortality (M1).
#' Returns an FLQuant with the specified constant value matching the dimensions
#' of the natural mortality from the input object.
#'
#' @param x An FLQuant object containing natural mortality, or an object with 
#'   a natural mortality method (e.g., FLStock, FLBRP)
#' @param y Numeric value for the background mortality rate (default = 0.025)
#' @return An FLQuant object with background mortality rates
#'
#' @details
#' M1 represents the constant or background component of natural mortality,
#' typically set to a small baseline value (default 0.025). The function creates
#' an FLQuant with this constant value matching the dimensions of the input.
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' library(FLBRP)
#' data(ple4brp)
#' # Using FLBRP object
#' M1=M1Fn(ple4brp)
#' M1_custom=M1Fn(ple4brp, y = 0.03)
#' 
#' # Using FLQuant object
#' m_vals=m(ple4brp)
#' M1=M1Fn(m_vals)
#' }
#'
#' @export
setGeneric("M1Fn", function(x, y = 0.025) standardGeneric("M1Fn"))

#' @rdname M1Fn
#' @export
setMethod("M1Fn", signature(x = "FLQuant"),
  function(x, y = 0.025) {
    FLQuant(y, dimnames = dimnames(x))
  })

#' @rdname M1Fn
#' @export
setMethod("M1Fn", signature(x = "FLStock"),
  function(x, y = 0.025) {
    FLQuant(y, dimnames = dimnames(m(x)))
  })

#' @rdname M1Fn
#' @export
setMethod("M1Fn", signature(x = "FLBRP"),
  function(x, y = 0.025) {
    FLQuant(y, dimnames = dimnames(m(x)))
  })

#' Calculate Predation/Variable Mortality (M2)
#'
#' @description
#' Calculates the variable/predation component of natural mortality (M2)
#' as the difference between total natural mortality (M) and background
#' mortality (M1).
#'
#' @param x An FLQuant object containing natural mortality, or an object with 
#'   natural mortality method (e.g., FLStock, FLBRP)
#' @param y Optional numeric value for M1 (default = 0.025)
#' @return An FLQuant object with M2 mortality rates (M - M1)
#'
#' @details
#' M2 represents the variable component of natural mortality, often associated
#' with predation or other density-dependent factors. It is calculated as:
#' M2 = M - M1
#' where M is total natural mortality and M1 is the background mortality.
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' library(FLBRP)
#' data(ple4brp)
#' # Using FLBRP object
#' M2=M2Fn(ple4brp)
#' 
#' # Using FLQuant object
#' m_vals=m(ple4brp)
#' M2=M2Fn(m_vals)
#' }
#'
#' @export
setGeneric("M2Fn", function(x, y = 0.025) standardGeneric("M2Fn"))

#' @rdname M2Fn
#' @export
setMethod("M2Fn", signature(x = "FLQuant"),
  function(x, y = 0.025) {
    M1=M1Fn(x, y)
    x - M1
  })

#' @rdname M2Fn
#' @export
setMethod("M2Fn", signature(x = "FLStock"),
  function(x, y = 0.025) {
    m(x) - M1Fn(x, y)
  })

#' @rdname M2Fn
#' @export
setMethod("M2Fn", signature(x = "FLBRP"),
  function(x, y = 0.025) {
    m(x) - M1Fn(x, y)
  })

#' Calculate Forage Biomass
#'
#' @description
#' Calculates forage biomass, representing the biomass available to predators,
#' weighted by mortality components and total mortality.
#'
#' @param x An FLQuant object or an object with stock data (e.g., FLStock, FLBRP)
#' @return An FLQuant object with forage biomass by year and iteration
#'
#' @details
#' Forage biomass is calculated as:
#' Forage = sum over ages of [stock.wt * stock.n * (M - M1) / Z * (1 - exp(-Z))]
#' where:
#' - M is natural mortality
#' - M1 is background mortality
#' - Z is total mortality (M + F)
#'
#' This represents the biomass available to predators, accounting for
#' the variable mortality component (M2 = M - M1) and mortality rates.
#'
#' For FLQuant objects, this method is not directly applicable as it requires
#' multiple stock components (stock.wt, stock.n, m, z). Use with FLStock or FLBRP objects instead.
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' library(FLBRP)
#' data(ple4brp)
#' forage_bio=forage(ple4brp)
#' }
#'
#' @export
setGeneric("forage", function(x) standardGeneric("forage"))

#' @rdname forage
#' @export
setMethod("forage", signature(x = "FLQuant"),
  function(x) {
    stop("forage method for FLQuant requires additional context (stock.wt, stock.n, m, z). Use forage(FLStock) or forage(FLBRP) instead.")
  })

#' @rdname forage
#' @export
setMethod("forage", signature(x = "FLStock"),
  function(x) {
    zVal=m(x) %+% harvest(x)
    temp = stock.wt(x) %*% stock.n(x) %*% (m(x) - M1Fn(x)) %/% zVal %*% (1 - exp(-zVal))
    res = FLCore::apply(temp, c(2, 6), sum)
    res
  })

#' @rdname forage
#' @export
setMethod("forage", signature(x = "FLBRP"),
  function(x) {
    zVal=m(x) %+% harvest(x)
    temp = stock.wt(x) %*% stock.n(x) %*% (m(x) %-% M1Fn(x)) %/% zVal %*%
            (1 - exp(-zVal))
    res = FLCore::apply(temp, c(2, 6), sum)
    res
  })
