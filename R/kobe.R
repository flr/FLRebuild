# =============================================================================
# Kobe Plot Functions
# =============================================================================

utils::globalVariables(c("bandwidth.nrd", "HPDregionplot", "mcmc"))

#' Calculate Kobe Quadrant Probabilities (HPD Region Method)
#'
#' @description
#' Calculates the probability of observations occurring in a 2D cell using
#' HPDregionplot. Given a sample, calculates the bivariate region of highest
#' marginal posterior density for two variables, using kde2d from MASS to
#' calculate a bivariate density.
#'
#' @param x A numeric vector
#' @param y A numeric vector
#' @param prob Numeric vector of probability levels (default = c(0.5, 0.75, 0.95))
#' @param n Number of points at which to evaluate the density grid (default = 21)
#' @param h Bandwidth of 2D kernel smoother. If not specified, defaults to
#'   bandwidth.nrd() values for x and y
#' @param lims Limits, specified as (x.lower, x.upper, y.lower, y.upper)
#'   (passed to kde2d)
#' @param na.rm Logical; if TRUE, any NA and NaN's are removed from x and y
#'   before calculations (default = FALSE)
#'
#' @return A data.frame with three variables:
#'   \item{x}{x coordinates of the grid points, vector of length n}
#'   \item{y}{y coordinates of the grid points, vector of length n}
#'   \item{level}{Contours corresponding to prob}
#'
#' @details
#' This function uses HPDregionplot from the coda package to calculate
#' highest posterior density regions for bivariate data. It's an alternative
#' implementation to the simpler probFn function which uses discrete thresholds.
#'
#' @examples
#' \dontrun{
#' x <- rnorm(20)
#' y <- rnorm(20)
#' probFn2(x, y)
#' }
#'
#' @keywords internal
probFn2 <- function(x, y, prob = c(0.5, 0.75, 0.95), n = 21,
                    h = c(bandwidth.nrd(x), bandwidth.nrd(y)),
                    lims = NULL, na.rm = FALSE) {
  
  if (na.rm) {
    .na <- !(is.na(x) | is.na(y) | is.nan(x) | is.nan(y))
    x <- x[.na]
    y <- y[.na]
  }
  
  tmp <- HPDregionplot(mcmc(data.frame(x, y)), prob = prob, h = h)
  
  prb <- ldply(tmp, function(dat) data.frame(level = dat$level, x = dat$x, y = dat$y))
  
  return(prb)
}

#' Calculate Kobe Quadrant Probabilities (Threshold Method)
#'
#' @description
#' Calculates probabilities for Kobe plot quadrants based on threshold values.
#' Determines whether stock and harvest are above/below reference points (1.0)
#' and calculates probabilities for red (collapsed), green (sustainable),
#' yellow (overfished only), and orange (overfishing only) quadrants.
#'
#' @param x A data.frame with columns 'stock' and 'harvest' containing
#'   relative values (e.g., B/BMSY, F/FMSY)
#'
#' @return A data.frame with columns:
#'   \item{red}{Probability of being in red quadrant (collapsed: B < 1 and F > 1)}
#'   \item{green}{Probability of being in green quadrant (sustainable: B >= 1 and F <= 1)}
#'   \item{yellow}{Probability of being in yellow quadrant (overfished only: B < 1 and F <= 1)}
#'   \item{orange}{Probability of being in orange quadrant (overfishing only: B >= 1 and F > 1)}
#'   \item{overFished}{Probability that stock is overfished (B < 1)}
#'   \item{overFishing}{Probability that overfishing is occurring (F > 1)}
#'
#' @details
#' The function converts stock and harvest values to binary indicators:
#' \itemize{
#'   \item b = 1 if stock >= 1 (not overfished), 0 otherwise
#'   \item f = 1 if harvest <= 1 (not overfishing), 0 otherwise
#' }
#'
#' Then calculates:
#' \itemize{
#'   \item green = f * b (both sustainable)
#'   \item red = (1 - b) * (1 - f) (both unsustainable)
#'   \item yellow = 1 - green - red (other combinations)
#' }
#'
#' @keywords internal
probFn <- function(x) {
  
  b <- pmax(pmin(as.integer(x$stock >= 1), 1), 0)
  f <- pmax(pmin(as.integer(x$harvest <= 1), 1), 0)
  p <- f * b
  collapsed <- (1 - b) * (1 - f)
  
  red   <- collapsed
  green <- p
  yellow <- 1 - p - collapsed
  
  overFished  <- 1 - b
  overFishing <- 1 - f
  
  data.frame(red = red, green = green, yellow = yellow,
             overFished = overFished, overFishing = overFishing)
}

#' Calculate Kobe Plot Probabilities
#'
#' @description
#' Calculates probabilities for Kobe plot quadrants from stock and harvest values.
#' Provides methods for calculating quadrant probabilities using threshold-based methods.
#'
#' @param stock Numeric vector or data.frame. If numeric, harvest must also be provided.
#'   If data.frame, should have columns 'stock' and 'harvest'.
#' @param harvest Numeric vector of harvest values (required if stock is numeric)
#'
#' @return A data.frame with quadrant probabilities (red, green, yellow, orange,
#'   overFished, overFishing)
#'
#' @details
#' This function provides methods for calculating Kobe quadrant probabilities
#' using threshold-based methods. The function determines which quadrant each
#' observation falls into based on whether stock and harvest are above/below
#' reference points (1.0).
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' # Using numeric vectors
#' stock_vals <- c(0.8, 1.2, 0.9, 1.1)
#' harvest_vals <- c(0.9, 0.8, 1.1, 1.2)
#' prob(stock_vals, harvest_vals)
#'
#' # Using data.frame
#' df <- data.frame(stock = stock_vals, harvest = harvest_vals)
#' prob(df)
#' }
#'
#' @export
setGeneric("prob", function(stock, harvest, ...) standardGeneric("prob"))

#' @rdname prob
#' @export
setMethod('prob', signature(stock = "numeric", harvest = "numeric"),
          function(stock, harvest) {
            
            res <- probFn(data.frame(stock = stock, harvest = harvest))
            
            return(res)
          })

#' @rdname prob
#' @export
setMethod('prob', signature(stock = 'data.frame', harvest = "missing"),
          function(stock) {
            
            # Ensure required columns exist
            if (!all(c("stock", "harvest") %in% names(stock))) {
              stop("data.frame must contain 'stock' and 'harvest' columns")
            }
            
            res <- probFn(stock)
            
            return(res)
          })

#' Calculate Kobe Plot Summary Statistics
#'
#' @description
#' Calculates summary statistics for Kobe plots from FLQuant objects or data.frames
#' containing stock biomass and harvest rate time series. Returns year-wise averages
#' of quadrant probabilities (red, green, yellow, orange).
#'
#' @param stock An FLQuant object or data.frame containing relative stock biomass values
#'   (e.g., B/BMSY). If data.frame, must contain 'stock' and 'harvest' columns.
#' @param harvest An FLQuant object containing relative harvest rate values (e.g., F/FMSY).
#'   Not required if stock is a data.frame.
#' @param ... Additional arguments
#'
#' @return A data.frame with columns:
#'   \item{year}{Year (if available)}
#'   \item{red}{Mean probability of being in red quadrant (collapsed)}
#'   \item{green}{Mean probability of being in green quadrant (sustainable)}
#'   \item{yellow}{Mean probability of being in yellow quadrant (overfished only)}
#'   \item{orange}{Mean probability of being in orange quadrant (overfishing only)}
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Combines stock and harvest FLQuants/data into a data.frame
#'   \item Calculates quadrant probabilities using prob() function
#'   \item Adds orange and yellow indicators (overfishing only vs overfished only)
#'   \item Aggregates to year-wise averages across iterations (if year column exists)
#' }
#'
#' Orange indicates overfishing but stock not overfished (B >= 1, F > 1).
#' Yellow indicates stock overfished but not overfishing (B < 1, F <= 1).
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' library(FLBRP)
#' data(ple4)
#' data(ple4brp)
#'
#' # Calculate relative stock and harvest
#' stock_rel <- ssb(ple4) %/% refpts(ple4brp)["msy", "ssb"]
#' harvest_rel <- fbar(ple4) %/% refpts(ple4brp)["msy", "harvest"]
#'
#' # Calculate Kobe summary from FLQuant objects
#' kobe_summary <- kobeSummary(stock_rel, harvest_rel)
#'
#' # Or from data.frame
#' df <- data.frame(stock = c(0.8, 1.2), harvest = c(0.9, 1.1), year = c(2010, 2011))
#' kobe_summary <- kobeSummary(df)
#' }
#'
#' @export
setGeneric("kobeSummary",
           function(stock, harvest, ...) standardGeneric("kobeSummary"))

#' @rdname kobeSummary
#' @export
setMethod("kobeSummary", signature(stock = "data.frame", harvest = "missing"),
  function(stock) {
    # Handle data.frame input directly
    if (!all(c("stock", "harvest") %in% names(stock))) {
      stop("data.frame must contain 'stock' and 'harvest' columns")
    }
    
    dtK <- stock
    
    ## Add quadrant probabilities
    dtK <- cbind(dtK, probFn(dtK))
    
    ## Orange: overfishing only; yellow: overfished only
    dtK <- dplyr::mutate(
      dtK,
      orange = as.numeric(overFishing == 1 & overFished == 0),
      yellow = as.numeric(overFishing == 0 & overFished == 1)
    )
    
    ## Year-wise averages (if year column exists)
    if ("year" %in% names(dtK)) {
      dtK <- plyr::ddply(
        dtK, .(year), with,
        data.frame(
          red    = mean(red),
          green  = mean(green),
          yellow = mean(yellow),
          orange = mean(orange)
        )
      )
    }
    
    dtK
  })

#' @rdname kobeSummary
#' @export
setMethod("kobeSummary",
          signature(stock = "FLQuant", harvest = "FLQuant"),
          function(stock, harvest) {
            
            ## Build FLQuants and data.frame
            dtK <- model.frame(FLQuants(stock = stock, harvest = harvest))
            
            ## Add quadrant probabilities (prob expects columns 'stock','harvest')
            dtK <- cbind(dtK, probFn(dtK))
            
            ## Orange: overfishing only; yellow: overfished only
            dtK <- dplyr::mutate(
              dtK,
              orange = as.numeric(overFishing == 1 & overFished == 0),
              yellow = as.numeric(overFishing == 0 & overFished == 1)
            )
            
            ## Year-wise averages
            dtK <- plyr::ddply(
              dtK, .(year), with,
              data.frame(
                red    = mean(red),
                green  = mean(green),
                yellow = mean(yellow),
                orange = mean(orange)
              )
            )
            
            dtK
          })
