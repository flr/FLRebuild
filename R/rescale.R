#' Rescale Numeric Values to a New Range
#'
#' @description
#' Linearly rescales numeric values from their current range to a specified
#' new range [a, b]. This is useful for normalizing data or converting values
#' to a specific scale.
#'
#' @param x Numeric vector to be rescaled
#' @param a Numeric. Lower bound of the target range (default range minimum)
#' @param b Numeric. Upper bound of the target range (default range maximum)
#'
#' @return A numeric vector of the same length as \code{x}, with values
#'   rescaled to the range [a, b]. If all values in \code{x} are identical
#'   (or if x_max == x_min), the function will return a vector of \code{a}
#'   values with a warning.
#'
#' @details
#' The rescaling is performed using the linear transformation:
#' \deqn{y = a + (x - x_{min}) \times \frac{b - a}{x_{max} - x_{min}}}{y = a + (x - x_min) * (b - a) / (x_max - x_min)}
#'
#' where \eqn{x_{min}}{x_min} and \eqn{x_{max}}{x_max} are the minimum and
#' maximum values in \code{x} (excluding NA values).
#'
#' Missing values (NA) are excluded from the min/max calculation but are
#' preserved in the output vector.
#'
#' @examples
#' # Rescale values to [0, 1] (normalization)
#' x <- c(10, 20, 30, 40, 50)
#' rescale(x, 0, 1)
#' # [1] 0.00 0.25 0.50 0.75 1.00
#'
#' # Rescale to [0, 100] (percentage scale)
#' rescale(x, 0, 100)
#' # [1]   0  25  50  75 100
#'
#' # Rescale to [-1, 1] (standardized range)
#' rescale(x, -1, 1)
#' # [1] -1.0 -0.5  0.0  0.5  1.0
#'
#' # Handle missing values
#' x_na <- c(10, 20, NA, 40, 50)
#' rescale(x_na, 0, 1)
#' # [1] 0.00 0.25   NA 0.75 1.00
#'
#' @seealso \code{\link{scale}} for standardization (z-scores)
#'
#' @export
rescale <- function(x, a, b) {
  # Validate inputs
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }
  if (!is.numeric(a) || length(a) != 1) {
    stop("'a' must be a single numeric value")
  }
  if (!is.numeric(b) || length(b) != 1) {
    stop("'b' must be a single numeric value")
  }
  if (a >= b) {
    stop("'a' must be less than 'b'")
  }
  
  # Calculate min and max (excluding NA)
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  
  # Handle case where all values are identical
  if (x_max == x_min) {
    warning("All values in 'x' are identical. Returning vector of 'a' values.")
    return(rep(a, length(x)))
  }
  
  # Perform linear rescaling
  result <- a + (x - x_min) * (b - a) / (x_max - x_min)
  
  return(result)
}
