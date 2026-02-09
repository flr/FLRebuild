#' Create Overlapping Time Windows from Time Series Data
#'
#' @description
#' Partitions time series data into overlapping windows based on year values.
#' This function is useful for rolling window analyses, where you want to
#' analyze data in overlapping time periods (e.g., 10-year windows with
#' 5-year overlap).
#'
#' @param dat A data.frame containing time series data
#' @param yearCol Character string. Name of the column containing year values
#'   (default: "year")
#' @param step Numeric. Step size between window centers in years (default: 10).
#'   Windows are centered at multiples of this step.
#' @param overlap Numeric. Number of years of overlap on each side of the
#'   window center (default: 5). The total window width is approximately
#'   \code{step + 2*overlap} years.
#'
#' @return A data.frame containing all windows stacked, with additional columns:
#'   \item{centre}{Numeric. The center year of each window}
#'   \item{from}{Numeric. The starting year of each window}
#'   \item{to}{Numeric. The ending year of each window}
#'   \item{.id}{Character. The center year as a character string (from plyr)}
#'
#' @details
#' The function creates windows as follows:
#' \enumerate{
#'   \item Determines the range of years in the data
#'   \item Calculates window centers at regular intervals (multiples of \code{step})
#'   \item For each center, creates a window from \code{centre - overlap} to
#'     \code{centre + step - 1 + overlap}
#'   \item Extracts all rows from \code{dat} that fall within each window
#'   \item Adds metadata columns (centre, from, to) to identify the window
#'   \item Combines all windows into a single data.frame
#' }
#'
#' Windows are overlapping, meaning the same data point can appear in multiple
#' windows. This is useful for:
#' \itemize{
#'   \item Rolling window regression analyses
#'   \item Time-varying parameter estimation
#'   \item Detecting regime shifts or trends over time
#'   \item Bootstrap resampling with temporal structure
#' }
#'
#' @examples
#' \dontrun{
#' # Create example time series data
#' dat <- data.frame(
#'   year = 1980:2020,
#'   value = rnorm(41),
#'   stock = rlnorm(41)
#' )
#'
#' # Create 10-year windows with 5-year overlap
#' windows <- windowOverlap(dat, yearCol = "year", step = 10, overlap = 5)
#'
#' # Analyze each window separately
#' library(dplyr)
#' windows %>%
#'   group_by(centre) %>%
#'   summarise(mean_value = mean(value), mean_stock = mean(stock))
#'
#' # Create 5-year windows with 2-year overlap
#' windows_small <- windowOverlap(dat, step = 5, overlap = 2)
#' }
#'
#' @seealso \code{\link[plyr]{ldply}} for list-to-data.frame conversion
#'
#' @importFrom plyr ldply
#' @export
windowOverlap <- function(dat, yearCol = "year", step = 10, overlap = 5) {
  # Validate inputs
  if (!is.data.frame(dat)) {
    stop("'dat' must be a data.frame")
  }
  if (!yearCol %in% names(dat)) {
    stop("Column '", yearCol, "' not found in 'dat'")
  }
  if (!is.numeric(step) || length(step) != 1 || step <= 0) {
    stop("'step' must be a positive numeric value")
  }
  if (!is.numeric(overlap) || length(overlap) != 1 || overlap < 0) {
    stop("'overlap' must be a non-negative numeric value")
  }
  
  # Extract years
  years <- dat[[yearCol]]
  
  if (!is.numeric(years)) {
    stop("Column '", yearCol, "' must contain numeric values")
  }
  
  # Calculate year range
  y_min <- min(years, na.rm = TRUE)
  y_max <- max(years, na.rm = TRUE)
  
  # Calculate window centers (multiples of step)
  centres <- seq(
    from = ceiling(y_min / step) * step,
    to   = floor(y_max / step) * step,
    by   = step
  )
  
  # If no valid centers, return empty data.frame with same structure
  if (length(centres) == 0) {
    warning("No valid window centers found. Returning empty data.frame.")
    result <- dat[FALSE, ]
    result$centre <- numeric(0)
    result$from <- numeric(0)
    result$to <- numeric(0)
    return(result)
  }
  
  # Create windows for each center
  out_list <- lapply(centres, function(y0) {
    # Define window boundaries
    from <- y0 - overlap
    to   <- y0 + step - 1 + overlap
    
    # Select rows within window
    sel <- years >= from & years <= to
    w   <- dat[sel, ]
    
    # Add window metadata
    w$centre <- y0
    w$from   <- from
    w$to     <- to
    
    return(w)
  })
  
  # Name the list by center year
  names(out_list) <- as.character(centres)
  
  # Combine all windows into single data.frame
  out <- plyr::ldply(out_list, .id = ".id")
  
  return(out)
}
