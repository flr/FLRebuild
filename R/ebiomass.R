#' Calculate Exploitable Biomass
#'
#' @description
#' Calculates the exploitable biomass from an FLStock object using selectivity-weighted
#' catch weights and stock numbers.
#'
#' @param object An object of class FLStock
#'
#' @return An FLQuant object containing the exploitable biomass time series
#'
#' @details
#' Exploitable biomass is based on weighting catch weights by selectivity normalised 
#' by peak selectivity, then multiplying by stock numbers and summing across ages
#'
#' @examples
#' \dontrun{
#' data(ple4)
#' eb <- ebiomass(ple4)
#' }
#'
#' @export
setGeneric("ebiomass", function(object) standardGeneric("ebiomass"))

#' Calculate total mortality-at-age (Z) using length data.
#' 
#' This is a generic S4 method that calculates something using the `haupt` function.
#' The specific implementation depends on the class of the input object.
#' 
#' @param object An object of class FLQuant or data.frame (depends on the method).
#' @param pars A parameter object (class FLPar) containing necessary parameters.
#' @param lc A threshold value for lc.
#' @param lmax A threshold value for lmax.
#' @param ... Additional arguments .
#' 
#' @return estimates of Z.
#' 
#' @export


#' @rdname ebiomass
#' @export
setMethod("ebiomass", signature(object="FLStock"),
          function(object) {
            sel   <- FLCore::harvest(object)
            wt    <- catch.wt(object) %*% sel %/% fapex(sel)
            eb.wt <- qmax(wt, 0.000001)
            
            apply(eb.wt %*% stock.n(object), 2:6, sum)
          })

setMethod("ebiomass", signature(object="FLBRP"),
          function(object) {
            sel   <- FLCore::harvest(object)
            wt    <- catch.wt(object) %*% sel %/% fapex(sel)
            eb.wt <- qmax(wt, 0.000001)
            
            apply(eb.wt %*% stock.n(object), 2:6, sum)
          })
