# =============================================================================
# Helper Functions
# =============================================================================

#' Interpolate rebuilding times
#' 
#' @description Finds the first year where biomass reaches the target (>= 1) for each initial depletion level
#' 
#' @param df Data frame with columns: initial, year, biomass
#' @return Data frame with columns: initial, year
#' @export
interp <- function(df) {
  # Input validation
  if (!is.data.frame(df) || !all(c("initial", "year", "biomass") %in% names(df))) {
    stop("df must be a data frame with columns: initial, year, biomass")
  }
  
  # Get unique initial values and sort them
  initials = sort(unique(df$initial))
  Yr = rep(NA, length(initials))
  
  # Find the first year where biomass >= 1 for each initial value
  for (i in seq_along(initials)) {
    sub_df = df[df$initial == initials[i], ]
    idx = which(sub_df$biomass >= 1)
    if (length(idx) > 0) {
      Yr[i] = sub_df$year[min(idx)]
    } else {
      Yr[i] = NA
    }
  }
  
  # Return result
  result = data.frame(
    initial = initials,
    year = Yr
  )
  
  return(result)
}

#' Safe try wrapper
#' 
#' @description Wraps a function call in try() and returns NULL on error instead of erroring
#' 
#' @param x Expression to try
#' @return Result of x if successful, NULL if error
#' @keywords internal
tryIt <- function(x) {
  rtn = try(x)
  if ("try-error" %in% is(rtn)) return(NULL)
  return(rtn)
}

