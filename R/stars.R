# =============================================================================
# Regime Shift Detection Using STARS Algorithm
# =============================================================================
# Based on Rodionov (2004) Sequential T-test Algorithm for Regime Shifts (STARS)
# Modified by Szuwalski et al. (submitted)
# =============================================================================

#' Detect Regime Shifts in Time Series Data
#'
#' @description
#' Detects regime shifts in a time series using a modified Sequential T-test
#' Algorithm for Regime Shifts (STARS; Rodionov 2004). This implementation
#' differs from the original STARS algorithm in that when a new regime is
#' detected, the algorithm begins searching at n+regimeLength (where n is the
#' first year of the regime identified), rather than n+1 as in the original.
#' Additionally, a 'Huber parameter' (which adjusts the influence of outliers)
#' is not included in this version.
#'
#' @param data Numeric vector of time series data to analyze
#' @param year Optional numeric vector of years corresponding to data points.
#'   If NULL (default), years are set to sequential indices (1, 2, 3, ...)
#' @param n Numeric. Minimum assumed regime length (default: 10). The algorithm
#'   requires at least this many observations to identify a regime.
#' @param sig Numeric. Significance level based on t-distribution (default: 1.68).
#'   This significance level is roughly equal to p < 0.1. Note: ideally this
#'   should be based on an appropriate t-distribution given the number of
#'   observations in the series, but this is not currently implemented.
#' @param plot Logical. If TRUE, plots the time series with regime polygons
#'   overlaid (default: FALSE)
#'
#' @return A data.frame with columns:
#'   \item{regime}{Factor identifying the regime number}
#'   \item{data}{Numeric values for plotting regime polygons (mean ± SD)}
#'   \item{year}{Numeric year values corresponding to regime boundaries}
#'
#' @details
#' The function iteratively identifies regime shifts in the time series:
#' \enumerate{
#'   \item Identifies the first regime shift using \code{ROregimeSHFT()}
#'   \item Continues searching for subsequent shifts starting from the end
#'     of the previous regime plus the minimum regime length
#'   \item Calculates mean and standard deviation for each identified regime
#'   \item Returns a data frame suitable for plotting regime polygons
#' }
#'
#' The returned data frame is structured for use with ggplot2, where each
#' regime is represented as a polygon with vertices at (minyear, mean+sd),
#' (minyear, mean-sd), (maxyear, mean-sd), (maxyear, mean+sd).
#'
#' @references
#' Rodionov, S. N. (2004). A sequential algorithm for testing climate regime
#' shifts. Geophysical Research Letters, 31(9), L09204.
#'
#' @examples
#' \dontrun{
#' # Create example time series with regime shifts
#' set.seed(123)
#' val <- c(rnorm(12, 0, 0.5), rnorm(12, 1, 0.5), rnorm(12, 0, 0.5))
#' years <- 1980:(1980 + length(val) - 1)
#'
#' # Detect regime shifts
#' regimes <- rodFn(val, year = years, n = 10, plot = TRUE)
#'
#' # Plot with ggplot2
#' library(ggplot2)
#' df <- data.frame(year = years, data = val)
#' ggplot(df, aes(x = year, y = data)) +
#'   geom_polygon(aes(year, data, group = regime),
#'                fill = "lavender", col = "blue",
#'                lwd = 0.25, data = regimes, alpha = 0.75) +
#'   geom_point() +
#'   geom_line()
#' }
#'
#' @seealso \code{\link{rod}} for FLQuant method, \code{\link{ROregimeSHFT}}
#'   for the internal shift detection algorithm
#'
#' @export
rodFn <- function(data, year = NULL, n = 10, sig = 1.68, plot = FALSE) {
  
  # Set default years if not provided
  if (is.null(year)) {
    year <- seq_along(data)
  }
  
  # Validate inputs
  if (length(data) != length(year)) {
    stop("Length of 'data' and 'year' must be equal")
  }
  
  if (length(data) <= n) {
    warning("Time series length (", length(data), 
            ") is less than or equal to minimum regime length (", n, 
            "). Cannot detect regime shifts.")
    return(data.frame(regime = factor(1),
                      data = c(mean(data) + sd(data), mean(data) - sd(data),
                               mean(data) - sd(data), mean(data) + sd(data)),
                      year = c(year[1], year[1], year[length(year)], 
                               year[length(year)])))
  }
  
  # Set the assumed minimum regime length
  regime_length <- n
  
  # Significance level (ideally should be based on t-distribution)
  significance <- sig
  
  # Iteration for finding all shifts in the time series
  shifts <- rep(NA, 100)  # Vector for recording shifts
  
  # Find first shift
  if (length(data) > regime_length) {
    shifts[1] <- tryCatch(
      ROregimeSHFT(regLN = regime_length, sig = significance, series = data),
      error = function(e) NA
    )
    
    # Continue finding subsequent shifts
    count <- 2
    if (!is.na(shifts[count - 1])) {
      while (shifts[count - 1] > 0 && 
             (shifts[count - 1] + regime_length + regime_length - 1) < length(data)) {
        shifts[count] <- ROregimeSHFT(
          regLN = regime_length,
          sig = significance,
          series = data,
          shift = shifts[count - 1] + regime_length - 1
        )
        count <- count + 1
      }
    }
  }
  
  # Extract valid shifts
  shifts_vec <- shifts[!is.na(shifts)]
  
  # Plot time series if requested
  if (plot) {
    plot(data, type = "l", xlab = "Index", ylab = "Value")
  }
  
  # Initialize vectors for regime statistics
  means <- NULL
  sds <- NULL
  lengths <- NULL
  
  # Calculate statistics for each regime
  for (i in seq_along(shifts_vec)) {
    if (i == 1) {
      # First regime: from start to first shift + regime_length - 1
      end_idx <- shifts_vec[i] + regime_length - 1
      reg_mean <- mean(data[1:end_idx])
      reg_sd <- sd(data[1:end_idx])
      means <- c(means, reg_mean)
      sds <- c(sds, reg_sd)
      lengths <- c(lengths, end_idx)
      
      if (plot) {
        polygon(
          x = c(1, end_idx, end_idx, 1),
          y = c(reg_mean + reg_sd, reg_mean + reg_sd,
                reg_mean - reg_sd, reg_mean - reg_sd),
          border = NA,
          col = "#0000ff55"
        )
      }
    } else {
      # Subsequent regimes: from previous shift end to current shift end
      start_idx <- shifts_vec[i - 1] + regime_length
      end_idx <- shifts_vec[i] + regime_length - 1
      
      # Handle last regime
      if (is.na(shifts_vec[i + 1])) {
        end_idx <- length(data)
      }
      
      reg_mean <- mean(data[start_idx:end_idx])
      reg_sd <- sd(data[start_idx:end_idx])
      means <- c(means, reg_mean)
      sds <- c(sds, reg_sd)
      lengths <- c(lengths, end_idx)
      
      if (plot) {
        polygon(
          x = c(start_idx, end_idx, end_idx, start_idx),
          y = c(reg_mean + reg_sd, reg_mean + reg_sd,
                reg_mean - reg_sd, reg_mean - reg_sd),
          border = NA,
          col = "#0000ff55"
        )
      }
    }
  }
  
  # If no shifts found, return single regime
  if (length(means) == 0) {
    means <- mean(data)
    sds <- sd(data)
    lengths <- length(data)
  }
  
  # Create regime data frame
  minyear <- c(year[1], year[rev(rev(lengths)[-1])] + 1)
  maxyear <- min(year) + lengths - 1
  
  regime_data <- data.frame(
    mn = means,
    sd = sds,
    i = factor(seq_along(means)),
    ln = lengths,
    minyear = minyear,
    maxyear = maxyear
  )
  
  # Format for plotting (polygon vertices)
  result <- data.frame(
    regime = regime_data$i,
    data = with(regime_data, c(mn + sd, mn - sd, mn - sd, mn + sd)),
    year = with(regime_data, c(minyear, minyear, maxyear, maxyear))
  )
  
  return(result)
}

#' Generic Function for Regime Shift Detection
#'
#' @description
#' Generic function for detecting regime shifts in time series data using the
#' STARS algorithm (Rodionov 2004) as modified by Szuwalski et al.
#'
#' @param object An object containing time series data (e.g., numeric vector, FLQuant)
#' @param ... Additional arguments passed to methods:
#'   \itemize{
#'     \item \code{year}: Optional year vector (defaults to sequential indices for numeric)
#'     \item \code{n}: Minimum regime length (default: 10)
#'     \item \code{sig}: Significance level (default: 1.68)
#'     \item \code{plot}: Logical, whether to plot (default: FALSE)
#'   }
#'
#' @return A data.frame with regime shift information
#'
#' @seealso \code{\link{rodFn}} for the base function, \code{\link{rod,numeric-method}}
#'   for numeric method, \code{\link{rod,FLQuant-method}} for FLQuant method
#'
#' @export
setGeneric("rod", function(object, year = NULL, n = 10, sig = 1.68, plot = FALSE, ...) 
           standardGeneric("rod"))

#' Detect Regime Shifts in Numeric Vectors
#'
#' @description
#' S4 method for detecting regime shifts in numeric vectors. This method
#' provides the same interface and defaults as the FLQuant method, ensuring
#' consistency across different input types.
#'
#' @param object A numeric vector containing time series data
#' @param year Optional numeric vector of years corresponding to data points.
#'   If NULL (default), years are set to sequential indices (1, 2, 3, ...)
#' @param n Numeric. Minimum assumed regime length (default: 10)
#' @param sig Numeric. Significance level based on t-distribution (default: 1.68)
#' @param plot Logical. If TRUE, plots the time series with regime polygons
#'   overlaid (default: FALSE)
#' @param ... Additional arguments (currently unused)
#'
#' @return A data.frame with regime shift information, structured for plotting
#'   with ggplot2. Contains columns: \code{regime}, \code{data}, \code{year}.
#'
#' @details
#' This method calls \code{rodFn} internally with the same default parameters
#' as the FLQuant method, ensuring consistent behavior across different input types.
#'
#' @examples
#' \dontrun{
#' # Create example time series with regime shifts
#' set.seed(123)
#' val <- c(rnorm(12, 0, 0.5), rnorm(12, 1, 0.5), rnorm(12, 0, 0.5))
#' years <- 1980:(1980 + length(val) - 1)
#'
#' # Detect regime shifts using rod method
#' regimes <- rod(val, year = years)
#'
#' # Plot with ggplot2
#' library(ggplot2)
#' df <- data.frame(year = years, data = val)
#' ggplot(df, aes(x = year, y = data)) +
#'   geom_polygon(aes(year, data, group = regime),
#'                fill = "lavender", col = "blue",
#'                lwd = 0.25, data = regimes, alpha = 0.75) +
#'   geom_point() +
#'   geom_line()
#' }
#'
#' @aliases rod,numeric-method
#' @rdname rod
#' @export
setMethod("rod", signature(object = "numeric"),
          function(object, year = NULL, n = 10, sig = 1.68, plot = FALSE, ...) {
            rodFn(data = object, year = year, n = n, sig = sig, plot = plot)
          })

#' Detect Regime Shifts in FLQuant Objects
#'
#' @description
#' S4 method for detecting regime shifts in \code{FLQuant} objects from the
#' FLCore package. This method applies \code{rodFn} to each iteration of the
#' FLQuant object separately, using the same default parameters as the numeric method.
#'
#' @param object An \code{FLQuant} object containing time series data
#' @param year Optional numeric vector of years. If NULL (default), years are
#'   extracted from the FLQuant object dimensions
#' @param n Numeric. Minimum assumed regime length (default: 10)
#' @param sig Numeric. Significance level based on t-distribution (default: 1.68)
#' @param plot Logical. If TRUE, plots the time series with regime polygons
#'   overlaid (default: FALSE)
#' @param ... Additional arguments (currently unused)
#'
#' @return A data.frame with regime shift information, structured for plotting
#'   with ggplot2. Contains columns: \code{regime}, \code{data}, \code{year},
#'   and \code{iter} (iteration number from FLQuant).
#'
#' @details
#' Evidence for regime shifts are explored using a sequential t-test algorithm
#' (STARS; Rodionov 2004) as modified by Szuwalski et al. (submitted).
#'
#' The function processes each iteration of the FLQuant object separately,
#' allowing for regime shift detection across multiple Monte Carlo iterations
#' or scenarios. The default parameters (n=10, sig=1.68, plot=FALSE) match those
#' of the numeric method for consistency.
#'
#' @references
#' Rodionov, S. N. (2004). A sequential algorithm for testing climate regime
#' shifts. Geophysical Research Letters, 31(9), L09204.
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' 
#' # Create example FLQuant with regime shifts
#' object <- FLQuant(rlnorm(1, FLQuant(0, dimnames = list(year = 1:30)), 0.3))
#' 
#' # Detect regime shifts (uses same defaults as numeric method)
#' regimes <- rod(object)
#' 
#' # Plot with ggplot2
#' library(ggplot2)
#' ggplot(as.data.frame(object)) +
#'   geom_polygon(aes(year, data, group = regime),
#'                fill = "lavender", col = "blue",
#'                lwd = 0.25, data = regimes, alpha = 0.75) +
#'   geom_point(aes(year, data)) +
#'   geom_line(aes(year, data))
#' }
#'
#' @aliases rod,FLQuant-method
#' @rdname rod
#' @export
setMethod("rod", signature(object = "FLQuant"),
          function(object, year = NULL, n = 10, sig = 1.68, plot = FALSE, ...) {
            plyr::ddply(
              as.data.frame(object),
              .(iter),
              function(df) {
                # Use years from FLQuant if year not provided
                if (is.null(year)) {
                  rodFn(data = df$data, year = df$year, n = n, sig = sig, plot = plot)
                } else {
                  rodFn(data = df$data, year = year, n = n, sig = sig, plot = plot)
                }
              }
            )
          })

#' Internal Function: Detect First Regime Shift in Time Series
#'
#' @description
#' Internal helper function that identifies the first regime shift in a time
#' series given a significance level and assumed regime length. This function
#' must be used iteratively (via \code{rodFn}) to identify all shifts in a
#' time series.
#'
#' @param regLN Numeric. Assumed minimum regime length (number of observations)
#' @param sig Numeric. Significance level based on t-distribution
#' @param series Numeric vector. Time series data to analyze
#' @param shift Numeric. Starting index for search (default: 0, meaning start
#'   from beginning). Used when searching for subsequent shifts.
#'
#' @return Numeric. The index (position) where the regime shift occurs, or 0
#'   if no shift is found. The returned value represents the length of the
#'   assumed regime before the shift occurs. For example, if the function
#'   returns 3 and the assumed regime length is 10, the regime shift occurs
#'   at index 13 (3 + 10).
#'
#' @details
#' This function implements a modified version of the STARS algorithm:
#' \itemize{
#'   \item When a new regime is detected, the algorithm begins searching at
#'     n + regimeLength (where n is the first year of the regime identified),
#'     whereas the original STARS begins searching at n + 1
#'   \item A 'Huber parameter' (which adjusts the influence of outliers) is
#'     not included in this version
#' }
#'
#' The algorithm works by:
#' \enumerate{
#'   \item Sliding a window of length \code{regLN} through the time series
#'   \item Calculating the mean and variance of each window
#'   \item Testing if the value at position \code{regLN + i} exceeds the
#'     threshold (basemean ± sig * sqrt(2 * basevar / regLN))
#'   \item If threshold exceeded, calculating a Regime Shift Index (RSI) and
#'     confirming the shift persists
#'   \item Returning the index of the confirmed shift, or 0 if none found
#' }
#'
#' @note
#' This is an internal function and is not exported. Users should call
#' \code{rodFn} or \code{rod} instead, which handle the iterative application
#' of this function.
#'
#' @references
#' Rodionov, S. N. (2004). A sequential algorithm for testing climate regime
#' shifts. Geophysical Research Letters, 31(9), L09204.
#'
#' @keywords internal
#' @noRd
ROregimeSHFT <- function(regLN, sig, series, shift = 0) {
  shifts <- 0
  start_year <- 1
  
  if (shift != 0) {
    start_year <- shift
  }
  
  # Iterate through possible shift positions
  for (i in seq(start_year, length(series) - regLN)) {
    # Define window for baseline statistics
    window <- series[i:(regLN + (i - 1))]
    base_mean <- mean(window, na.rm = TRUE)
    base_var <- var(window, na.rm = TRUE)
    threshold <- sig * sqrt(2 * base_var / regLN)
    
    # Test if value at regLN + i exceeds threshold
    # Note: This differs from original STARS which tests at i + 1
    test_value <- series[regLN + i]
    
    if (test_value > base_mean + threshold || test_value < base_mean - threshold) {
      # Calculate Regime Shift Index (RSI)
      rsi <- 0
      if (test_value > base_mean + threshold) {
        new_mean <- base_mean + threshold
        rsi <- (test_value - new_mean) / (regLN * sqrt(base_var))
        direction_down <- 0
      } else {
        new_mean <- base_mean - threshold
        rsi <- (new_mean - test_value) / (regLN * sqrt(base_var))
        direction_down <- 1
      }
      
      original_sign <- sign(rsi)
      counter <- min(regLN - 1, length(series) - (regLN + (i - 1)))
      rsi <- 0
      
      # Confirm shift persists
      for (j in seq_len(counter)) {
        shifts <- 0
        if (direction_down == 0) {
          rsi <- rsi + (series[regLN + i + j - 1] - new_mean) / 
                       (regLN * sqrt(base_var))
        } else {
          rsi <- rsi + (new_mean - series[regLN + i + j - 1]) / 
                       (regLN * sqrt(base_var))
        }
        
        # Break if RSI changes sign or becomes infinite
        if (sign(rsi) != original_sign || is.infinite(rsi)) {
          break
        }
        
        # If we've checked all positions, confirm the shift
        if (j == counter) {
          shifts <- i
          break
        }
      }
    }
    
    # If shift confirmed, stop searching
    if (shifts != 0) {
      break
    }
  }
  
  return(shifts)
}
