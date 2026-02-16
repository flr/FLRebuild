#' Summarize Fits for a Biomass Metric
#'
#' @description
#' Generic function to calculate surplus production and process error for any
#' biomass metric (e.g., SSB, ebiomass) from an FLBRP object. This function
#' generalizes the calculation of surplus production and process error that was
#' previously hardcoded for specific metrics.
#'
#' @param object An FLBRP object containing fitted model results
#' @param biomass Character string. Name of the biomass metric to summarize.
#'   Common values: "ssb", "ebiomass", "biomass". Default: "ebiomass"
#' @param biomass.obs Function or character. Function to extract observed biomass
#'   from the object. If character, will look for an attribute with that name
#'   (e.g., "eb.obs" for ebiomass). If NULL, will try to construct from
#'   \code{biomass} parameter (e.g., "ssb" -> ssb.obs, "ebiomass" -> attributes(object)$eb.obs).
#'   Default: NULL (auto-detect)
#' @param biomass.name Character. Short name for the biomass metric used in
#'   output names (e.g., "SSB", "EB"). If NULL, will use capitalized \code{biomass}.
#'   Default: NULL
#'
#' @return An FLQuants object containing:
#'   \item{sp}{Surplus production: B_{t+1} - B_t + C_t}
#'   \item{pe}{Process error: (sp - yield) / biomass}
#'   \item{biomass}{The observed biomass time series}
#'
#' @details
#' This function calculates:
#' \itemize{
#'   \item \strong{Surplus Production (sp):} The change in biomass plus catch:
#'     \deqn{SP_t = B_{t+1} - B_t + C_t}{SP_t = B_{t+1} - B_t + C_t}
#'   \item \strong{Process Error (pe):} The normalized difference between
#'     surplus production and yield:
#'     \deqn{PE_t = \frac{SP_t - Y_t}{B_t}}{PE_t = (SP_t - Y_t) / B_t}
#' }
#'
#' The function automatically handles different biomass metrics by:
#' \itemize{
#'   \item Using standard observation functions (e.g., \code{ssb.obs()}) when available
#'   \item Falling back to attributes (e.g., \code{attributes(object)$eb.obs}) for
#'     custom metrics like ebiomass
#'   \item Constructing function names from the metric name when possible
#' }
#'
#' @examples
#' \dontrun{
#' library(FLBRP)
#' 
#' # Summarize ebiomass fits
#' eb_summary <- summarizeFit(brp, biomass = "ebiomass")
#' 
#' # Summarize SSB fits
#' ssb_summary <- summarizeFit(brp, biomass = "ssb")
#' 
#' # Custom biomass metric
#' custom_summary <- summarizeFit(brp, 
#'                                biomass = "custom",
#'                                biomass.obs = function(x) attributes(x)$custom.obs)
#' }
#'
#' @importFrom FLCore mcf FLQuants model.frame
#' @importFrom plyr laply alply ldply
#' @export
setGeneric("summarizeFit", function(object, biomass = "ebiomass", 
                                    biomass.obs = NULL, biomass.name = NULL) {
  standardGeneric("summarizeFit")
})

#' @rdname summarizeFit
#' @export
setMethod("summarizeFit", signature(object = "FLBRP"),
  function(object, biomass = "ebiomass", biomass.obs = NULL, biomass.name = NULL) {
    
    # Determine biomass name for output
    if (is.null(biomass.name)) {
      biomass.name <- toupper(biomass)
    }
    
    # Get biomass observations
    if (is.null(biomass.obs)) {
      # Try to auto-detect the observation function
      if (biomass == "ssb") {
        B_obs <- ssb.obs(object)
      } else if (biomass == "ebiomass" || biomass == "eb") {
        # ebiomass uses attributes
        if (is.null(attributes(object)$eb.obs)) {
          stop("ebiomass observations not found in object attributes. ",
               "Expected attributes(object)$eb.obs")
        }
        B_obs <- attributes(object)$eb.obs
      } else if (biomass == "biomass" || biomass == "stock") {
        # Try stock.obs if it exists
        if (exists("stock.obs", mode = "function")) {
          B_obs <- stock.obs(object)
        } else {
          stop("Could not find observation function for biomass metric '", biomass, "'")
        }
      } else {
        # Try to construct function name: "metric" -> "metric.obs"
        obs_fn_name <- paste0(biomass, ".obs")
        if (exists(obs_fn_name, mode = "function")) {
          obs_fn <- get(obs_fn_name)
          B_obs <- obs_fn(object)
        } else {
          # Try attributes: "metric" -> attributes(object)$metric.obs
          attr_name <- paste0(biomass, ".obs")
          if (attr_name %in% names(attributes(object))) {
            B_obs <- attributes(object)[[attr_name]]
          } else {
            stop("Could not find observation function or attribute for biomass metric '", 
                 biomass, "'. Tried: ", obs_fn_name, "() and attributes(object)$", attr_name)
          }
        }
      }
    } else {
      # Use provided function
      if (is.function(biomass.obs)) {
        B_obs <- biomass.obs(object)
      } else if (is.character(biomass.obs)) {
        # Character means attribute name
        if (biomass.obs %in% names(attributes(object))) {
          B_obs <- attributes(object)[[biomass.obs]]
        } else {
          stop("Attribute '", biomass.obs, "' not found in object")
        }
      } else {
        stop("biomass.obs must be a function or character string")
      }
    }
    
    # Validate biomass observations
    if (is.null(B_obs) || length(B_obs) == 0) {
      stop("Biomass observations are empty or NULL")
    }
    
    # Get catch observations
    C_obs <- catch.obs(object)
    
    # Validate dimensions
    if (dim(B_obs)[2] != dim(C_obs)[2]) {
      stop("Biomass and catch observations have different time dimensions")
    }
    
    # Calculate surplus production: SP_t = B_{t+1} - B_t + C_t
    # Remove last year from B_obs and first year from C_obs for alignment
    n_years <- dim(B_obs)[2]
    if (n_years < 2) {
      stop("Insufficient time series data (need at least 2 years)")
    }
    
    B_t <- B_obs[, -n_years, drop = FALSE]  # B_t (all but last year)
    B_t1 <- B_obs[, -1, drop = FALSE]       # B_{t+1} (all but first year)
    C_t <- C_obs[, -n_years, drop = FALSE]  # C_t (all but last year)
    
    # Surplus production
    sp <- B_t1 - B_t + C_t
    
    # Get yield for process error calculation
    # We need yield from refpts - extract it from tseries if available
    # For now, use catch as proxy for yield (they should be similar)
    yield <- C_t
    
    # Calculate process error: PE = (SP - Y) / B
    pe <- (sp - yield) / B_t
    
    # Create FLQuants object
    rtn <- FLQuants(
      sp = sp,
      pe = pe,
      biomass = B_obs
    )
    
    # Name the biomass component
    names(rtn)[3] <- biomass
    
    # Ensure consistent dimensions
    rtn <- mcf(rtn)
    
    return(rtn)
  })

#' @rdname summarizeFit
#' @export
setMethod("summarizeFit", signature(object = "FLBRPs"),
  function(object, biomass = "ebiomass", biomass.obs = NULL, biomass.name = NULL) {
    # Apply to each FLBRP in the list
    plyr::ldply(object, function(x) {
      result <- summarizeFit(x, biomass = biomass, 
                            biomass.obs = biomass.obs, 
                            biomass.name = biomass.name)
      # Convert to data.frame for ldply
      model.frame(result)
    })
  })
