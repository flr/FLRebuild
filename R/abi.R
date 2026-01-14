# =============================================================================
# Age-Based Indicators (ABI)
# =============================================================================

#' Calculate Reference Age for Age-Based Indicator
#'
#' @description
#' Calculates the reference age at which a specified proportion (p) of the
#' cumulative stock numbers is reached at a given reference point (e.g., MSY).
#' This age is used as a threshold for age-based indicator calculations.
#'
#' @param object An FLBRP object
#' @param ref Character string specifying the reference point (default = "msy")
#' @param p Numeric value between 0 and 1 specifying the cumulative proportion
#'   threshold (default = 0.9)
#'
#' @return An FLQuant object containing reference ages by year and iteration
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Sets fishing mortality to the specified reference point (e.g., FMSY)
#'   \item Calculates cumulative proportions of stock numbers by age
#'   \item Finds the minimum age at which cumulative proportion exceeds p
#'   \item Returns the reference age + 1 (or maximum age if threshold not reached)
#' }
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' library(FLBRP)
#' data(ple4brp)
#' 
#' # Calculate reference age at FMSY with 90% threshold
#' ref_age=abiAge(ple4brp, ref = "msy", p = 0.9)
#' 
#' # Calculate reference age at F0.1 with 95% threshold
#' ref_age=abiAge(ple4brp, ref = "f0.1", p = 0.95)
#' }
#'
#' @rdname abiAge
#' @export
setMethod("abiAge",
          signature(object = "FLBRP"),
          function(object, ref = "msy", p = 0.9) {
            
            fbar(object)=as.FLQuant(
              refpts(object)[ref, "harvest", drop = TRUE],
              dimnames = list(iter = seq(dim(object)[6]))
            )
            
            stk.n=stock.n(object)[-1]
            
            cumN =FLCore::apply(stk.n, c(2, 6), cumsum) %/% quantSums(stk.n)
            agesQ=ages(stk.n)
            
            agesQ[cumN <= p]=NA
            
            FLCore::apply(
              agesQ, c(2:6),
              function(x) min(c(x + 1, dims(object)$max), na.rm = TRUE)
            )
          })

#' Calculate Proportion of Stock Numbers Above Reference Age at MSY
#'
#' @description
#' Calculates the proportion of stock numbers (P(N)) at or above the reference age
#' at a specified reference point (e.g., MSY). This represents the proportion of
#' older fish in the population at equilibrium.
#'
#' @param object An FLBRP object
#' @param ref Character string specifying the reference point (default = "msy")
#' @param p Numeric value between 0 and 1 specifying the cumulative proportion
#'   threshold for reference age calculation (default = 0.9)
#'
#' @return An FLQuant object containing proportions (0-1) by year and iteration
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Calculates the reference age using \code{abiAge(object, ref, p)}
#'   \item Determines which ages are at or above the reference age
#'   \item Calculates the proportion of stock numbers in these age classes
#'   \item Returns the proportion as an FLQuant
#' }
#'
#' This indicator reflects the age structure of the population at equilibrium
#' conditions (e.g., at MSY).
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' library(FLBRP)
#' data(ple4brp)
#' 
#' # Calculate P(N) at MSY with 90% threshold
#' pmsy=abiMsy(ple4brp, ref = "msy", p = 0.9)
#' }
#'
#' @rdname abiMsy
#' @export
setMethod("abiMsy",
          signature(object = "FLBRP"),
          function(object, ref = "msy", p = 0.9) {
            
            fbar(object)=as.FLQuant(
              refpts(object)[ref, "harvest", drop = TRUE],
              dimnames = list(iter = seq(dim(object)[6]))
            )
            
            A    =abiAge(object, ref, p)
            stk.n=stock.n(object)[-1]
            
            flag=FLQuant(
              ages(stk.n) >= FLCore::expand(A, age = dimnames(stk.n)$age)
            )
            
            FLCore::apply(stk.n %*% flag, c(2, 6), sum) %/% apply(stk.n, c(2, 6), sum)
          })

#' Calculate Proportion of Stock Numbers Above Threshold Age
#'
#' @description
#' Helper function that calculates the proportion of stock numbers at or above
#' a specified threshold age for an FLStock object.
#'
#' @param x An FLStock object
#' @param A An FLQuant object containing threshold age(s) by year and iteration
#'
#' @return An FLQuant object containing proportions (0-1) by year and iteration
#'
#' @details
#' This is an internal helper function used by the \code{abi} methods. It:
#' \enumerate{
#'   \item Extracts stock numbers (excluding age 0)
#'   \item Propagates threshold ages if dimensions don't match
#'   \item Creates a flag matrix for ages >= threshold age
#'   \item Calculates proportion of stock numbers above threshold
#' }
#'
#' @keywords internal
abiStock <- function(x, A) {
  
  stk.n=stock.n(x)[-1]
  
  if (dim(x)[6] > 1 && dim(A)[6] == 1)
    A=propagate(A, dim(x)[6])
  
  amsy=FLQuant(
    rep(c(A), each = prod(dim(stk.n)[-6])),
    dimnames = dimnames(stk.n)
  )
  
  flag=FLQuant(ages(stk.n) >= amsy)
  
  FLCore::apply(stk.n %*% flag, c(2, 6), sum) %/% apply(stk.n, c(2, 6), sum)
}

#' Calculate Age-Based Indicator (ABI)
#'
#' @description
#' Calculates the Age-Based Indicator, which compares the proportion of older fish
#' in the observed stock to the proportion expected at a reference point (e.g., MSY).
#' ABI values > 1 indicate a higher proportion of older fish than expected at the
#' reference point, while values < 1 indicate fewer older fish.
#'
#' @param object An FLStock object
#' @param age An FLBRP object (for reference point calculation) or FLQuant object
#'   (for direct threshold age specification)
#' @param ref Character string specifying the reference point when age is FLBRP
#'   (default = "msy")
#' @param p Numeric value between 0 and 1 specifying the cumulative proportion
#'   threshold when age is FLBRP (default = 0.9)
#'
#' @return An FLQuant object containing ABI values by year and iteration
#'
#' @details
#' When \code{age} is an FLBRP object:
#' \enumerate{
#'   \item Calculates reference proportion at MSY using \code{abiMsy(age, ref, p)}
#'   \item Calculates reference age using \code{abiAge(age, ref, p)}
#'   \item Calculates observed proportion above reference age using \code{abiStock(object, age)}
#'   \item Returns the ratio: observed proportion / reference proportion
#' }
#'
#' When \code{age} is an FLQuant object:
#' \enumerate{
#'   \item Uses the FLQuant directly as threshold ages
#'   \item Calculates observed proportion above threshold ages
#'   \item Returns the proportion (not normalized by reference)
#' }
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' library(FLBRP)
#' data(ple4)
#' data(ple4brp)
#' 
#' # Calculate ABI using FLBRP reference
#' abi_vals=abi(ple4, ple4brp, ref = "msy", p = 0.9)
#' 
#' # Calculate ABI using direct threshold age
#' ref_age=abiAge(ple4brp, ref = "msy", p = 0.9)
#' abi_vals=abi(ple4, ref_age)
#' }
#'
#' @rdname abi
#' @export
setMethod("abi",
          signature(object = "FLStock", age = "FLBRP"),
          function(object, age, ref = "msy", p = 0.9) {
            
            pmsy=abiMsy(age, ref, p)
            ageQ=abiAge(age, ref, p)
            
            pt=abiStock(object, ageQ)
            
            pt %/% pmsy
          })

#' @rdname abi
#' @export
setMethod("abi",
          signature(object = "FLStock", age = "FLQuant"),
          function(object, age) {
            
            stk.n=stock.n(object)[-1]
            
            if (dim(object)[6] > 1 & dim(age)[6] == 1)
              age=propagate(age, dim(object)[6])
            
            amsy=FLQuant(
              rep(c(age), each = prod(dim(stk.n)[-6])),
              dimnames = dimnames(stk.n)
            )
            
            flag=FLQuant(ages(stk.n) >= amsy)
            
            apply(stk.n %*% flag, c(2, 6), sum) %/% apply(stk.n, c(2, 6), sum)
          })

# =============================================================================
# Backwards-Compatible Wrapper Functions
# =============================================================================

#' Age-Based Indicator (Backwards-Compatible Wrapper)
#'
#' @description
#' Backwards-compatible wrapper for \code{abi()} function. Maintains the old
#' function name for code compatibility.
#'
#' @param x An FLStock object
#' @param y An FLBRP object (default: FLBRP(x) if x can be converted)
#' @param ref Character string specifying the reference point (default = "msy")
#' @param p Numeric value between 0 and 1 (default = 0.9)
#'
#' @return An FLQuant object containing ABI values
#'
#' @seealso \code{\link{abi}}
#'
#' @export
ABI <- function(x, y = FLBRP(x), ref = "msy", p = 0.9) {
  abi(x, y, ref = ref, p = p)
}

#' Age-Based Indicator at MSY (Backwards-Compatible Wrapper)
#'
#' @description
#' Backwards-compatible wrapper for \code{abiMsy()} function.
#'
#' @param y An FLBRP object
#' @param ref Character string specifying the reference point (default = "msy")
#' @param p Numeric value between 0 and 1 (default = 0.9)
#'
#' @return An FLQuant object containing proportions at MSY
#'
#' @seealso \code{\link{abiMsy}}
#'
#' @export
ABIMSY <- function(y, ref = "msy", p = 0.9) {
  abiMsy(y, ref = ref, p = p)
}

#' Age-Based Indicator for Stock (Backwards-Compatible Wrapper)
#'
#' @description
#' Backwards-compatible wrapper for \code{abi(x, A)} when A is an FLQuant.
#'
#' @param x An FLStock object
#' @param A An FLQuant object containing threshold age(s)
#'
#' @return An FLQuant object containing proportions
#'
#' @seealso \code{\link{abi}}
#'
#' @export
ABIstock <- function(x, A) {
  abi(x, A)
}

# =============================================================================
# ABI at Target F
# =============================================================================

#' Calculate Age-Based Indicator at Target Fishing Mortality
#'
#' @description
#' Calculates the Age-Based Indicator for a stock at a specified target fishing
#' mortality rate (Ftarget). The function creates an equilibrium stock at Ftarget,
#' calculates a reference proportion, and compares observed stock proportions to
#' this reference.
#'
#' @param stock An FLStock object
#' @param ftgt Numeric value specifying the target fishing mortality rate
#'   (default = 0.2)
#' @param thresh Numeric value between 0 and 1 specifying the cumulative proportion
#'   threshold for reference age calculation (default = 0.9)
#'
#' @return An FLQuant object containing ABI values normalized by the reference
#'   proportion at Ftarget
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Creates an FLBRP object from the stock
#'   \item Sets fishing mortality to near-zero for initialization, then to ftgt
#'   \item Calculates equilibrium stock at Ftarget
#'   \item Determines the reference age where cumulative proportion exceeds thresh
#'   \item Calculates reference proportion of stock numbers above reference age
#'   \item Compares observed stock proportions to reference and returns normalized ratio
#' }
#'
#' This provides an ABI metric relative to a specific target F rather than MSY.
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' data(ple4)
#' 
#' # Calculate ABI at F = 0.2
#' abi_tgt=ABItgt(ple4, ftgt = 0.2, thresh = 0.9)
#' 
#' # Calculate ABI at F = 0.3 with 95% threshold
#' abi_tgt=ABItgt(ple4, ftgt = 0.3, thresh = 0.95)
#' }
#'
#' @export
ABItgt <- function(stock, ftgt = 0.2, thresh = 0.9) {
  
  eqstk=brp(FLBRP(stock))
  
  fbar(eqstk)[, 1][]    =0.00001  # near-zero F
  fbar(eqstk)[, 1:101][]=ftgt     # equilibrium at ftgt
  
  eqstk=brp(eqstk)
  eqstk=FLCore::window(as(eqstk, "FLStock"), start = 2, end = 2)
  eqstk@name=stock@name
  
  n_a =stock.n(eqstk)[-1, ]
  ages=dims(n_a)$min:dims(n_a)$max
  
  cums    =FLCore::apply(n_a, 2:6, cumsum)
  n_thresh=sum(n_a * thresh)
  
  aref=min(
    ages[which((n_thresh - cums)^2 ==
                 min((n_thresh - cums)^2))] + 1,
    range(eqstk)["plusgroup"] - 1
  )
  
  rp=sum(n_a[ac(aref:max(ages)), ]) / sum(n_a)
  
  flq=quantSums(stock.n(stock)[ac(aref:range(stock)[2]), ]) /
    quantSums(stock.n(stock)[-1, ]) / rp
  
  flq
}
