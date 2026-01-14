#' Convert Logit-Transformed Steepness to Original Scale
#'
#' Converts steepness parameter from logit-transformed scale back to the original
#' scale used in fisheries models.
#'
#' @param logitH Numeric vector of logit-transformed steepness values
#' @return Numeric vector of steepness values on original scale (0.2 to 1.0)
#' 
#' @examples
#' # Convert logit steepness values to original scale (0.2 to 1.0)
#' logit_values <- c(0, 1, 2)
#' steepness <- fromLogits(logit_values)
#' print(steepness)
#' 
#' @export
fromLogits <- function(logitH) {
  # Validate input
  if (!is.numeric(logitH)) {
    stop("logitH must be numeric")
  }
  
  # Convert using inverse logit transformation with bounds
  out <- 0.2001 + 0.7998 * 1 / (1 + exp(-logitH))
  
  return(out)
}

#' Convert Steepness to Logit Scale
#'
#' Converts steepness parameter to logit-transformed scale for numerical
#' stability in optimization routines.
#'
#' @param h Numeric vector of steepness values (0.2 to 1.0)
#' @return Numeric vector of logit-transformed steepness values
#' 
#' @examples
#' # Convert steepness values to logit scale for optimization
#' steepness_values <- c(0.3, 0.5, 0.8)
#' logit_values <- toLogits(steepness_values)
#' print(logit_values)
#' 
#' @export
toLogits <- function(h) {
  # Validate input
  if (!is.numeric(h)) {
    stop("h must be numeric")
  }
  
  # Check bounds
  if (any(h <= 0.2001 | h >= 1.0)) {
    stop("h must be between 0.2001 and 1.0")
  }
  
  # Convert to logit scale
  result <- -log(0.7998 / (h - 0.2001) - 1)
  
  return(result)
}

