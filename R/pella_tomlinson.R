# Note: FLCore should be in package dependencies
# No need to load here - package imports handle this

#' Pella-Tomlinson Production Function S4 Class
#' 
#' @description
#' Implementation of the Pella-Tomlinson production function for fisheries modeling.
#' The production function is: P(B) = (r/p) * B * [1 - (B/virgin)^p]
#' 
#' @details
#' This package provides S4 classes and methods for working with Pella-Tomlinson
#' production models, including parameter estimation, reference point calculations,
#' and utility functions for data manipulation and plotting.
#' 
#' @author Your Name
#' @docType package
#' @name pellaTomlinson
NULL

#' Pella-Tomlinson Production Model Class
#' 
#' @description
#' S4 class representing a Pella-Tomlinson production model with parameters
#' and stock biomass data.
#' 
#' @slot params FLPar object containing model parameters (r, p, virgin, bmsy, msy)
#' @slot stock FLQuant object containing stock biomass data
#' 
#' @export
setClass("PellaTomlinson",
         slots = list(
           params = "FLPar",      # Model parameters (r, p, virgin, bmsy, msy)
           stock = "FLQuant"      # Stock biomass data
         ),
         validity = function(object) {
           # Check that required parameters exist
           required_params = c("r", "p", "virgin", "bmsy", "msy")
           param_names = unlist(dimnames(object@params)$params)
           if (!all(required_params %in% param_names)) {
             return("Missing required parameters: r, p, virgin, bmsy, msy")
           }
           
           # Check that virgin biomass is positive
           if (any(object@params["virgin"] <= 0)) {
             return("Virgin biomass must be positive")
           }
           
           # Check that BMSY/virgin ratio is valid (0 < BMSY/virgin < 1)
           ratio = object@params["bmsy"] %/% object@params["virgin"]
           if (any(ratio <= 0) || any(ratio >= 1)) {
             return("BMSY/virgin ratio must be between 0 and 1")
           }
           
           return(TRUE)
         })

#' Constructor for Pella-Tomlinson Objects
#' 
#' @description
#' Creates a new Pella-Tomlinson production model object.
#' 
#' @param params FLPar object containing model parameters
#' @param stock FLQuant object containing stock biomass data (optional)
#' 
#' @return A Pella-Tomlinson object
#' 
#' @details
#' The constructor automatically handles parameter estimation:
#' - If params contains bmsy, k (virgin), and msy, it estimates r and p
#' - If params contains r, p, and virgin, it uses them directly
#' - The k parameter is automatically renamed to virgin for consistency
#' 
#' @examples
#' # Create with MSY, BMSY, and K - automatically estimates r and p
#' params = FLPar(bmsy = 0.4, k = 1, msy = 0.1)
#' pt_model = PellaTomlinson(params)
#' 
#' # Create with r, p, and virgin (direct usage)
#' params = FLPar(r = 0.5, p = 0.25, virgin = 1.0, bmsy = 0.7, msy = 0.1)
#' pt_model = PellaTomlinson(params)
#' 
#' # Create with custom stock data
#' custom_stock = FLQuant(seq(0, 2, length.out = 50))
#' pt_model = PellaTomlinson(params, stock = custom_stock)
#' 
#' @export
PellaTomlinson = function(params, stock = NULL, ...) {
  # Handle case where individual parameters are passed as named arguments
  if (missing(params) && length(list(...)) > 0) {
    # Extract named arguments
    args = list(...)
    if (all(c("msy", "bmsy", "virgin") %in% names(args))) {
      # Create FLPar from individual parameters
      params = FLPar(
        msy = args$msy,
        bmsy = args$bmsy,
        virgin = args$virgin
      )
    } else if (all(c("r", "p", "virgin") %in% names(args))) {
      # Create FLPar from individual parameters
      params = FLPar(
        r = args$r,
        p = args$p,
        virgin = args$virgin
      )
    } else {
      stop("Must provide either msy, bmsy, virgin OR r, p, virgin as named arguments")
    }
  }
  
  # Check if we need to estimate r and p
  if (inherits(params, "FLPar")) {
    param_names = dimnames(params)$params
    
    # Case 1: If we have bmsy, k (or virgin), and msy, estimate r and p
    if (all(c("bmsy", "msy") %in% param_names) && 
        any(c("k", "virgin") %in% param_names)) {
      # Create a copy with k renamed to virgin for consistency
      params_renamed = params
      if ("k" %in% param_names) {
        pnames = dimnames(params_renamed)$params
        pnames[pnames == "k"] = "virgin"
        dimnames(params_renamed)$params = pnames
      }
      
                    # Handle multi-iteration FLPar objects
       dims = dim(params_renamed)
       n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
       

       
       if (n_iter == 1) {
        # Single iteration - use existing logic
        params_numeric = setNames(c(params_renamed), dimnames(params_renamed)$params)
        params_complete = param(params_numeric)
        params_complete_flpar = FLPar(params_complete)
        
        # Create default stock if not provided
        if (is.null(stock)) {
          virgin_val = c(params_complete_flpar["virgin"])
          stock = FLQuant(seq(0, virgin_val, length.out = 101))
        }
        
        return(new("PellaTomlinson", params = params_complete_flpar, stock = stock))
      } else {
        # Multiple iterations - estimate r and p for each
        bmsy_vals = c(params_renamed["bmsy"])
        virgin_vals = c(params_renamed["virgin"])
        msy_vals = c(params_renamed["msy"])
        
        r_vals = numeric(n_iter)
        p_vals = numeric(n_iter)
        
        for (i in 1:n_iter) {
          params_i = c(bmsy = bmsy_vals[i], virgin = virgin_vals[i], msy = msy_vals[i])
          p_i = pFn(params_i)
          r_i = rFn(params_i)
          p_vals[i] = p_i
          r_vals[i] = r_i
        }
        
        # Create matrix-based FLPar for multiple iterations
        param_matrix = matrix(
          c(bmsy_vals, virgin_vals, msy_vals, p_vals, r_vals),
          nrow = 5,
          ncol = n_iter,
          byrow = TRUE,
          dimnames = list(
            params = c("bmsy", "virgin", "msy", "p", "r"),
            iter = 1:n_iter
          )
        )
        params_complete = FLPar(param_matrix)
        
        # Create default stock for multi-iteration case
        if (is.null(stock)) {
          # Create stock with max virgin value to cover all iterations
          max_virgin = max(virgin_vals)
          stock = FLQuant(seq(0, max_virgin, length.out = 101))
        }
        
        return(new("PellaTomlinson", params = params_complete, stock = stock))
      }
    }
    # Case 2: If we already have r, p, and virgin, compute bmsy & msy and use directly
    else if (all(c("r", "p", "virgin") %in% param_names)) {
      # Handle multi-iteration FLPar objects
      dims = dim(params)
      n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
      
      if (n_iter == 1) {
        # Single iteration - use existing logic
        # Ensure params has bmsy and msy; if missing, derive and append
        pnames = dimnames(params)$params
        need_bmsy = !("bmsy" %in% pnames)
        need_msy  = !("msy" %in% pnames)
        if (need_bmsy || need_msy) {
          # Work with numeric named vector to compute
          pn = setNames(c(params), pnames)
          if (need_bmsy) pn["bmsy"] = bmsy(pn)
          if (need_msy)  pn["msy"]  = msy(pn)
          params = FLPar(pn)
        }
        
        # Create default stock if not provided
        if (is.null(stock)) {
          virgin_val = c(params["virgin"])
          stock = FLQuant(seq(0, virgin_val, length.out = 101))
        }
        
        return(new("PellaTomlinson", params = params, stock = stock))
      } else {
        # Multiple iterations - compute bmsy and msy for each
        r_vals = c(params["r"])
        p_vals = c(params["p"])
        virgin_vals = c(params["virgin"])
        
        bmsy_vals = numeric(n_iter)
        msy_vals = numeric(n_iter)
        
        for (i in 1:n_iter) {
          params_i = c(r = r_vals[i], p = p_vals[i], virgin = virgin_vals[i])
          bmsy_vals[i] = bmsy(params_i)
          msy_vals[i] = msy(params_i)
        }
        
                 if (n_iter == 1) {
           # Single iteration - create simple FLPar
           params_complete = FLPar(
             r = r_vals[1],
             p = p_vals[1],
             virgin = virgin_vals[1],
             bmsy = bmsy_vals[1],
             msy = msy_vals[1]
           )
         } else {
           # Multiple iterations - create matrix-based FLPar
           param_matrix = matrix(
             c(r_vals, p_vals, virgin_vals, bmsy_vals, msy_vals),
             nrow = n_iter,
             ncol = 5,
             byrow = TRUE,
             dimnames = list(
               iter = 1:n_iter,
               params = c("r", "p", "virgin", "bmsy", "msy")
             )
           )
           params_complete = FLPar(param_matrix)
         }
        
        # Create default stock for multi-iteration case
        if (is.null(stock)) {
          # Create stock with max virgin value to cover all iterations
          max_virgin = max(virgin_vals)
          stock = FLQuant(seq(0, max_virgin, length.out = 101))
        }
        
        return(new("PellaTomlinson", params = params_complete, stock = stock))
      }
    }
         # Case 3: If we have bmsy, virgin, and msy, estimate r and p
     else if (all(c("bmsy", "virgin", "msy") %in% param_names)) {
       # Handle multi-iteration FLPar objects
       dims = dim(params)
       n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
       
       if (n_iter == 1) {
         # Single iteration - use existing logic
         params_numeric = setNames(c(params), dimnames(params)$params)
         params_complete = param(params_numeric)
         params_complete_flpar = FLPar(params_complete)
         
         # Create default stock if not provided
         if (is.null(stock)) {
           virgin_val = c(params_complete_flpar["virgin"])
           stock = FLQuant(seq(0, virgin_val, length.out = 101))
         }
         
         return(new("PellaTomlinson", params = params_complete_flpar, stock = stock))
       } else {
         # Multiple iterations - estimate r and p for each
         bmsy_vals = c(params["bmsy"])
         virgin_vals = c(params["virgin"])
         msy_vals = c(params["msy"])
         
         r_vals = numeric(n_iter)
         p_vals = numeric(n_iter)
         
         for (i in 1:n_iter) {
           params_i = c(bmsy = bmsy_vals[i], virgin = virgin_vals[i], msy = msy_vals[i])
           p_i = pFn(params_i)
           r_i = rFn(params_i)
           p_vals[i] = p_i
           r_vals[i] = r_i
         }
         
         # Create complete FLPar with all iterations
         # We need to create a properly structured FLPar with iterations
         # Create a matrix with parameters as rows and iterations as columns
         param_matrix = matrix(
           c(bmsy_vals, virgin_vals, msy_vals, p_vals, r_vals),
           nrow = 5,
           ncol = n_iter,
           byrow = TRUE,
           dimnames = list(
             params = c("bmsy", "virgin", "msy", "p", "r"),
             iter = 1:n_iter
           )
         )
         
         # Create FLPar with proper iteration structure
         params_complete = FLPar(param_matrix)
         
         # Create default stock for multi-iteration case
         if (is.null(stock)) {
           # Create stock with max virgin value to cover all iterations
           max_virgin = max(virgin_vals)
           stock = FLQuant(seq(0, max_virgin, length.out = 101))
         }
         
         return(new("PellaTomlinson", params = params_complete, stock = stock))
       }
     }
         # Case 4: Invalid parameter combination
     else {
       stop("params must contain either:\n",
            "  - bmsy, k (or virgin), and msy (will estimate r and p)\n",
            "  - r, p, and virgin (direct usage)")
     }
  } else {
    stop("params must be an FLPar object")
  }
}

#' Accessor method for parameters
#' 
#' @description
#' Retrieves the parameters from a Pella-Tomlinson object.
#' 
#' @param object Pella-Tomlinson object
#' @param ... Additional arguments
#' 
#' @return FLPar object containing the model parameters
#' 
#' @examples
#' params = FLPar(r = 0.5, p = 0.25, virgin = 1.0, bmsy = 0.7, msy = 0.1)
#' pt_model = PellaTomlinson(params)
#' model_params = params(pt_model)
#' 
#' @export
setGeneric("params", function(object, ...) standardGeneric("params"))

#' @rdname params
#' @export
setMethod("params", signature(object = "PellaTomlinson"), 
          function(object) {
            return(object@params)
          })

setMethod("params", signature(object = "numeric"), 
          function(object) {
            # For numeric vectors, return a named numeric vector
            if (length(object) == 0) return(numeric(0))
            
            # Try to get names from the object
            param_names = names(object)
            if (is.null(param_names)) {
              # If no names, try to infer from length and common patterns
              if (length(object) == 3) {
                param_names = c("r", "p", "virgin")
              } else if (length(object) == 5) {
                param_names = c("r", "p", "virgin", "bmsy", "msy")
              } else {
                param_names = paste0("param", seq_along(object))
              }
            }
            
            return(setNames(object, param_names))
          })

#' Initialize Parameter Vector
#' 
#' @description
#' Creates a parameter vector with the specified values for Pella-Tomlinson models.
#' 
#' @param msy Maximum Sustainable Yield (optional)
#' @param bmsy Biomass at MSY (optional)
#' @param virgin Virgin biomass/carrying capacity (default: 1)
#' @param r Intrinsic growth rate (optional)
#' @param p Shape parameter (optional)
#' 
#' @return Named numeric vector with parameters
#' 
#' @examples
#' # Initialize with MSY and BMSY
#' params = initParams(msy = 0.1, bmsy = 0.7, virgin = 1.0)
#' 
#' # Initialize with r and p
#' params = initParams(r = 0.5, p = 0.25, virgin = 1.0)
#' 
#' @export
initParams = function(msy = NA, bmsy = NA, virgin = 1, r = NA, p = NA) {
  c(r = r, p = p, virgin = virgin, bmsy = bmsy, msy = msy)
}

#' Calculate Shape Parameter p from BMSY and Virgin Biomass
#' 
#' @description
#' Numerically solves for the shape parameter p given BMSY and virgin biomass.
#' Uses the relationship: BMSY = virgin * (1/(1+p))^(1/p)
#' 
#' @param params Named numeric vector containing 'bmsy' and 'virgin' parameters
#' 
#' @return Numeric value of the shape parameter p
#' 
#' @details
#' The function uses uniroot to solve the non-linear equation for p.
#' For p near zero, it uses the Fox model approximation (p = 0).
#' 
#' @examples
#' params = c(bmsy = 0.7, virgin = 1.0)
#' p_val = pFn(params)
#' 
#' @export
pFn = function(params) {
  bmsy = params["bmsy"]
  virgin = params["virgin"]

  # Check if Bmsy/virgin ratio is valid
  ratio = bmsy / virgin
  if (ratio <= 0 || ratio >= 1) {
    warning("Bmsy/virgin ratio must be between 0 and 1")
    return(NA)
  }
  
  # Use numerical approach to find p
  fn = function(p, shape) {
    if (abs(p) < 1e-10) {
      return(exp(-1) - shape)
    } else {
      return((1 / (1 + p))^(1 / p) - shape)
    }
  }
  
  # Try to find p using uniroot
  tryCatch({
    result = uniroot(fn, shape = ratio, lower = -0.999, upper = 100, 
                     tol = 1e-8, extendInt = "upX")
    return(result$root)
  }, error = function(e) {
    warning("Numerical solution failed, returning default p = 0.25")
    return(0.25)
  })
}

#' Calculate Intrinsic Growth Rate r from MSY, BMSY, and Virgin Biomass
#' 
#' @description
#' Calculates the intrinsic growth rate r given MSY, BMSY, and virgin biomass.
#' 
#' @param params Named numeric vector containing 'msy', 'bmsy', and 'virgin' parameters
#' 
#' @return Numeric value of the intrinsic growth rate r
#' 
#' @details
#' First calculates p using pFn(), then solves for r using the relationship:
#' MSY = (r/p) * BMSY * (1 - (BMSY/virgin)^p)
#' 
#' @examples
#' params = c(msy = 0.1, bmsy = 0.7, virgin = 1.0)
#' r_val = rFn(params)
#' 
#' @export
rFn = function(params) {
  msy = params["msy"]
  bmsy = params["bmsy"]
  virgin = params["virgin"]
  
  # First get p from Bmsy
  p_val = pFn(params)
  if (is.na(p_val)) return(NA)
  
  # For the Pella-Tomlinson model, solve for r directly
  # MSY = (r/p) * BMSY * (1 - (BMSY/virgin)^p)
  # Rearranging: r = MSY * p / (BMSY * (1 - (BMSY/virgin)^p))
  r_val = msy * p_val / (bmsy * (1 - (bmsy / virgin)^p_val))
  
  return(r_val)
}

#' Estimate Both r and p Parameters
#' 
#' @description
#' Estimates both r and p parameters from MSY, BMSY, and virgin biomass.
#' 
#' @param params Named numeric vector containing 'msy', 'bmsy', and 'virgin' parameters
#' 
#' @return Updated parameter vector with estimated 'r' and 'p' values
#' 
#' @examples
#' params = c(msy = 0.1, bmsy = 0.7, virgin = 1.0)
#' estimated_params = param(params)
#' 
#' @export
param = function(params) {
  # Update p in params
  p_val = pFn(params)
  if (is.na(p_val)) {
    params["p"] = NA
    params["r"] = NA
    return(params)
  }
  
  params["p"] = p_val
  
  # Update r in params
  r_val = rFn(params)
  params["r"] = r_val
  
  return(params)
}

#' Production Function
#' 
#' @description
#' Calculates the production (yield) given biomass and model parameters.
#' 
#' @param biomass Numeric vector of biomass values.
#' @param params Named numeric vector containing 'r', 'p', and 'virgin' parameters.
#' 
#' @return Numeric vector of production values.
#' 
#' @examples
#' biomass = seq(0, 1, length.out = 101)
#' params = c(r = 0.5, p = 0.25, virgin = 1.0)
#' yield = production(biomass, params)
#' 
#' @export
production = function(biomass, params) {
  r = params["r"]
  p = params["p"]
  virgin = params["virgin"]
  
  # Handle negative p values carefully
  if (p < 0) {
    ratio = biomass / virgin
    if (any(ratio <= 0)) {
      warning("Biomass must be positive for negative p values")
      return(rep(NA, length(biomass)))
    }
  }
  
  (r / p) * biomass * (1 - (biomass / virgin)^p)
}

#' Calculate Maximum Sustainable Yield (MSY)
#' 
#' @description
#' Calculates the Maximum Sustainable Yield (MSY) given model parameters.
#' 
#' @param params Named numeric vector containing 'r', 'p', and 'virgin' parameters.
#' 
#' @return Numeric value of MSY.
#' 
#' @examples
#' params = c(r = 0.5, p = 0.25, virgin = 1.0)
#' msy_val = msy(params)
#' 
#' @export
msy = function(params) {
  r = params["r"]
  p = params["p"]
  virgin = params["virgin"]
  # For Pella-Tomlinson: MSY = (r/p) * BMSY * (1 - (BMSY/virgin)^p)
  bmsy = virgin * (1 / (1 + p))^(1 / p)
  (r / p) * bmsy * (1 - (bmsy / virgin)^p)
}

#' Calculate Biomass at MSY (BMSY)
#' 
#' @description
#' Calculates the Biomass at MSY (BMSY) given model parameters.
#' 
#' @param params Named numeric vector containing 'p' and 'virgin' parameters.
#' 
#' @return Numeric value of BMSY.
#' 
#' @examples
#' params = c(p = 0.25, virgin = 1.0)
#' bmsy_val = bmsy(params)
#' 
#' @export
bmsy = function(params) {
  p = params["p"]
  virgin = params["virgin"]
  # For Pella-Tomlinson: BMSY = virgin * (1/(1+p))^(1/p)
  virgin * (1 / (1 + p))^(1 / p)
}

#' Calculate Fishing Mortality at MSY (FMSY)
#' 
#' @description
#' Calculates the Fishing Mortality at MSY (FMSY) given model parameters.
#' 
#' @param params Named numeric vector containing 'r' and 'p' parameters.
#' 
#' @return Numeric value of FMSY.
#' 
#' @examples
#' params = c(r = 0.5, p = 0.25)
#' fmsy_val = fmsy(params)
#' 
#' @export
fmsy = function(params) {
  r = params["r"]
  p = params["p"]
  r / (p+1)
}

#' S4 Method Definitions for the PellaTomlinson class
#' 
#' @description
#' Generic functions for the PellaTomlinson class.
#' 
#' @details
#' These functions are used to define the methods for the S4 class.
#' 
#' @name PellaTomlinson-methods
#' @keywords internal
NULL

#' Generic function for production
#' 
#' @description
#' Generic function for calculating production.
#' 
#' @param biomass Numeric vector or FLQuant object.
#' @param params FLPar object.
#' @param ... Additional arguments.
#' 
#' @return Numeric vector or FLQuant object.
#' 
#' @export
setGeneric("production", function(biomass, params, ...) standardGeneric("production"))

#' Generic function for MSY
#' 
#' @description
#' Generic function for calculating Maximum Sustainable Yield (MSY).
#' 
#' @param params FLPar object.
#' @param ... Additional arguments.
#' 
#' @return Numeric value or FLQuant object.
#' 
#' @export
setGeneric("msy", function(params, ...) standardGeneric("msy"))

#' Generic function for BMSY
#' 
#' @description
#' Generic function for calculating Biomass at MSY (BMSY).
#' 
#' @param params FLPar object.
#' @param ... Additional arguments.
#' 
#' @return Numeric value or FLQuant object.
#' 
#' @export
setGeneric("bmsy", function(params, ...) standardGeneric("bmsy"))

#' Generic function for FMSY
#' 
#' @description
#' Generic function for calculating Fishing Mortality at MSY (FMSY).
#' 
#' @param params FLPar object.
#' @param ... Additional arguments.
#' 
#' @return Numeric value or FLQuant object.
#' 
#' @export
setGeneric("fmsy", function(params, ...) standardGeneric("fmsy"))

#' Generic function for p
#' 
#' @description
#' Generic function for calculating the shape parameter p.
#' 
#' @param params FLPar object.
#' @param ... Additional arguments.
#' 
#' @return FLPar object with p parameter.
#' 
#' @export
setGeneric("p", function(params, ...) standardGeneric("p"))

#' Generic function for r
#' 
#' @description
#' Generic function for calculating the intrinsic growth rate r.
#' 
#' @param params FLPar object.
#' @param ... Additional arguments.
#' 
#' @return FLPar object with r parameter.
#' 
#' @export
setGeneric("r", function(params, ...) standardGeneric("r"))

#' Generic function for param
#' 
#' @description
#' Generic function for estimating both r and p parameters.
#' 
#' @param params FLPar object.
#' @param ... Additional arguments.
#' 
#' @return FLPar object with estimated parameters.
#' 
#' @export
setGeneric("param", function(params, ...) standardGeneric("param"))





#' Methods for numeric params
#' 
#' @description
#' Methods for calculating production, MSY, BMSY, FMSY, p, r, param
#' when parameters are provided as numeric vectors.
#' 
#' @name numeric-methods
#' @keywords internal
NULL

#' Method for production (numeric)
#' 
#' @description
#' Method for calculating production when parameters are provided as a numeric vector.
#' 
#' @param biomass Numeric vector.
#' @param params Numeric vector.
#' 
#' @return Numeric vector.
#' 
#' @export
setMethod("production", signature(biomass = "numeric", params = "numeric"), 
          function(biomass, params) {
            r = params["r"]
            p = params["p"]
            virgin = params["virgin"]
            
            # Handle negative p values carefully
            if (p < 0) {
              ratio = biomass / virgin
              if (any(ratio <= 0)) {
                warning("Biomass must be positive for negative p values")
                return(rep(NA, length(biomass)))
              }
            }
            
            (r / p) * biomass * (1 - (biomass / virgin)^p)
          })

#' Method for MSY (numeric)
#' 
#' @description
#' Method for calculating MSY when parameters are provided as a numeric vector.
#' 
#' @param params Numeric vector.
#' 
#' @return Numeric value.
#' 
#' @export
setMethod("msy", signature(params = "numeric"), 
          function(params) {
            r = params["r"]
            p = params["p"]
            virgin = params["virgin"]
            # For Pella-Tomlinson: MSY = (r/p) * BMSY * (1 - (BMSY/virgin)^p)
            bmsy = virgin * (1 / (1 + p))^(1 / p)
            (r / p) * bmsy * (1 - (bmsy / virgin)^p)
          })

#' Method for BMSY (numeric)
#' 
#' @description
#' Method for calculating BMSY when parameters are provided as a numeric vector.
#' 
#' @param params Numeric vector.
#' 
#' @return Numeric value.
#' 
#' @export
setMethod("bmsy", signature(params = "numeric"), 
          function(params) {
            p = params["p"]
            virgin = params["virgin"]
            # For Pella-Tomlinson: BMSY = virgin * (1/(1+p))^(1/p)
            virgin * (1 / (1 + p))^(1 / p)
          })

#' Method for FMSY (numeric)
#' 
#' @description
#' Method for calculating FMSY when parameters are provided as a numeric vector.
#' 
#' @param params Numeric vector.
#' 
#' @return Numeric value.
#' 
#' @export
setMethod("fmsy", signature(params = "numeric"), 
          function(params) {
            r = params["r"]
            p = params["p"]
            r / p
          })

#' Method for p (numeric)
#'
#' @description
#' Method for calculating the shape parameter p when parameters are provided as a numeric vector.
#'
#' @param params Numeric vector.
#'
#' @return Numeric value of the shape parameter p.
#'
#' @export
setMethod("p", signature(params = "numeric"),
          function(params) {
            # Use the robust pFn utility function for consistency
            p_val = pFn(params)
            return(p_val)
          })

#' Method for r (numeric)
#' 
#' @description
#' Method for calculating the intrinsic growth rate r when parameters are provided as a numeric vector.
#' 
#' @param params Numeric vector.
#' 
#' @return Numeric value of the intrinsic growth rate r.
#' 
#' @export
setMethod("r", signature(params = "numeric"), 
          function(params) {
            msy = params["msy"]
            bmsy = params["bmsy"]
            virgin = params["virgin"]
            
            p_val = pFn(params)
            if (is.na(p_val)) return(NA_real_)
            
            # For the Pella-Tomlinson model, solve for r directly
            r_val = msy * p_val / (bmsy * (1 - (bmsy / virgin)^p_val))
            
            return(r_val)
          })

#' Method for param (numeric)
#' 
#' @description
#' Method for estimating both r and p parameters when parameters are provided as a numeric vector.
#' 
#' @param params Numeric vector.
#' 
#' @return Numeric vector with estimated parameters.
#' 
#' @export
setMethod("param", signature(params = "numeric"), 
          function(params) {
            p_val = pFn(params)
            if (is.na(p_val)) {
              params["p"] = NA
              params["r"] = NA
              return(params)
            }
            
            params["p"] = p_val
            r_val = rFn(params)
            params["r"] = r_val
            
            return(params)
          })

#' Methods for FLPar objects
#' 
#' @description
#' Methods for calculating production, MSY, BMSY, FMSY, p, r, param
#' when parameters are provided as FLPar objects.
#' 
#' @name FLPar-methods
#' @keywords internal
NULL

#' Method for production (FLPar)
#' 
#' @description
#' Method for calculating production when parameters are provided as an FLPar object.
#' 
#' @param biomass FLQuant object.
#' @param params FLPar object.
#' 
#' @return FLQuant object.
#' 
#' @export
setMethod("production", signature(biomass = "FLQuant", params = "FLPar"), 
          function(biomass, params) {
            r = params["r"]
            p = params["p"]
            virgin = params["virgin"]
            ((r %/% p) %*% biomass) %*% (1 - exp(log(biomass %/% virgin) %*% p))
          })

#' Method for MSY (FLPar)
#' 
#' @description
#' Method for calculating MSY when parameters are provided as an FLPar object.
#' 
#' @param params FLPar object.
#' 
#' @return Numeric value or FLQuant object.
#' 
#' @export
setMethod("msy", signature(params = "FLPar"), 
          function(params) {
            r = params["r"]
            p = params["p"]
            virgin = params["virgin"]
            # For Pella-Tomlinson: MSY = (r/p) * BMSY * (1 - (BMSY/virgin)^p)
            bmsy = virgin %*% exp(log(1 / (1 + p)) %*% (1 / p))
            res = ((r %/% p) %*% bmsy) %*% (1 - exp(log(bmsy %/% virgin) %*% p))
            if (inherits(res, "FLPar")) {
              dn = dimnames(res)
              dn$params[1] = "msy"
              dimnames(res) = dn
              return(res)
            } else {
              return(FLPar(msy = c(res)))
            }
          })

#' Method for BMSY (FLPar)
#' 
#' @description
#' Method for calculating BMSY when parameters are provided as an FLPar object.
#' 
#' @param params FLPar object.
#' 
#' @return Numeric value or FLQuant object.
#' 
#' @export
setMethod("bmsy", signature(params = "FLPar"), 
          function(params) {
            p = params["p"]
            virgin = params["virgin"]
            # For Pella-Tomlinson: BMSY = virgin * (1/(1+p))^(1/p)
            res = virgin %*% exp(log(1 / (1 + p)) %*% (1 / p))
            if (inherits(res, "FLPar")) {
              dn = dimnames(res)
              dn$params[1] = "bmsy"
              dimnames(res) = dn
              return(res)
            } else {
              return(FLPar(bmsy = c(res)))
            }
          })

#' Method for FMSY (FLPar)
#' 
#' @description
#' Method for calculating FMSY when parameters are provided as an FLPar object.
#' 
#' @param params FLPar object.
#' 
#' @return Numeric value or FLQuant object.
#' 
#' @export
setMethod("fmsy", signature(params = "FLPar"), 
          function(params) {
            r = params["r"]
            p = params["p"]
            r %/% p
          })

#' Method for p (FLPar)
#'
#' @description
#' Method for calculating the shape parameter p when parameters are provided as an FLPar object.
#'
#' @param params FLPar object.
#'
#' @return FLPar object with p parameter.
#'
#' @export
setMethod("p", signature(params = "FLPar"),
          function(params) {
            # Use the robust pFn utility function to avoid duplication and ensure consistency
            bmsy = c(params["bmsy"])
            virgin = c(params["virgin"])
            
            # Create a named vector for pFn
            params_numeric = c(bmsy = bmsy, virgin = virgin)
            p_val = pFn(params_numeric)
            
            return(FLPar(p = p_val))
          })

#' Method for r (FLPar)
#' 
#' @description
#' Method for calculating the intrinsic growth rate r when parameters are provided as an FLPar object.
#' 
#' @param params FLPar object.
#' 
#' @return FLPar object with r parameter.
#' 
#' @export
setMethod("r", signature(params = "FLPar"), 
          function(params) {
            msy = c(params["msy"])
            bmsy = c(params["bmsy"])
            virgin = c(params["virgin"])
            
            # Use pFn utility function to avoid circular dependency
            p_val = pFn(c(bmsy = bmsy, virgin = virgin))
            if (is.na(p_val)) return(FLPar(r = NA))
            
            # For the Pella-Tomlinson model, solve for r directly
            r_val = (msy * p_val) / (bmsy * (1 - (bmsy/virgin)^p_val))
            
            return(FLPar(r = r_val))
          })









#' Method for param (FLPar)
#'
#' @description
#' Method for estimating both r and p parameters when parameters are provided as an FLPar object.
#'
#' @param params FLPar object.
#'
#' @return FLPar object with estimated parameters.
#'
#' @export
setMethod("param", signature(params = "FLPar"),
          function(params) {
            # Use utility functions to avoid circular dependency
            p_val = pFn(c(bmsy = c(params["bmsy"]), virgin = c(params["virgin"])))
            if (is.na(p_val)) {
              # Create new FLPar with NA values for p and r
              return(FLPar(
                msy = c(params["msy"]),
                bmsy = c(params["bmsy"]),
                virgin = c(params["virgin"]),
                p = NA,
                r = NA
              ))
            }
            
            r_val = rFn(c(msy = c(params["msy"]), bmsy = c(params["bmsy"]), virgin = c(params["virgin"])))
            
            # Create new FLPar with all parameters
            return(FLPar(
              msy = c(params["msy"]),
              bmsy = c(params["bmsy"]),
              virgin = c(params["virgin"]),
              p = p_val,
              r = r_val
            ))
          })

#' Method for param (PellaTomlinson)
#' 
#' @description
#' Method for estimating both r and p parameters when parameters are provided as a PellaTomlinson object.
#' 
#' @param params PellaTomlinson object.
#' 
#' @return FLPar object with estimated parameters.
#' 
#' @export
setMethod("param", signature(params = "PellaTomlinson"), 
          function(params) {
            # Extract the FLPar from the PellaTomlinson object and call the FLPar method
            return(param(params@params))
          })

#' Methods for PellaTomlinson class
#' 
#' @description
#' Methods for calculating production, MSY, BMSY, FMSY, p, r, param
#' when the PellaTomlinson object is used as the method's first argument.
#' 
#' @name PellaTomlinson-methods
#' @keywords internal
NULL

#' Method for production (PellaTomlinson)
#' 
#' @description
#' Method for calculating production when the PellaTomlinson object is used.
#' 
#' @param biomass Missing (object's stock data is used).
#' @param params PellaTomlinson object.
#' 
#' @return FLQuant object.
#' 
#' @export
setMethod("production", signature(biomass = "PellaTomlinson", params = "missing"), 
          function(biomass, params) {
            # Use the stock data from the object
            biomass_data = biomass@stock
            params_obj = biomass@params
            
            # Check if we have multiple iterations
            dims = dim(params_obj)
            n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
            
            if (n_iter == 1) {
              # Single iteration - use existing logic
              r = c(params_obj["r"])
              p = c(params_obj["p"])
              virgin = c(params_obj["virgin"])
              return(((r / p) * biomass_data) * (1 - exp(log(biomass_data / virgin) * p)))
            } else {
              # Multiple iterations - handle each iteration properly
              # Create a list to store results for each iteration
              result_list = list()
              
              # Extract all parameters for all iterations
              r_vals = c(params_obj["r"])
              p_vals = c(params_obj["p"])
              virgin_vals = c(params_obj["virgin"])
              
              for (i in 1:n_iter) {
                # Extract parameters for this iteration
                r_i = r_vals[i]
                p_i = p_vals[i]
                virgin_i = virgin_vals[i]
                
                # Calculate production for this iteration
                prod_i = ((r_i / p_i) * biomass_data) * (1 - exp(log(biomass_data / virgin_i) * p_i))
                
                # Store result with iteration dimension
                result_list[[i]] = prod_i
              }
              
              # Combine results into a single FLQuant with iterations
              # Build an FLQuant explicitly over the iter dimension to avoid dplyr::combine masking
              proto = result_list[[1]]
              dims_proto = dim(proto)
              dn_proto = dimnames(proto)
              dims_proto[6] = length(result_list)
              dn_proto$iter = as.character(seq_len(length(result_list)))
              result = FLQuant(array(NA, dim = dims_proto, dimnames = dn_proto))
              for (i in seq_along(result_list)) {
                result[, , , , , i] = result_list[[i]]
              }
              
              return(result)
            }
          })



#' Method for production (PellaTomlinson)
#' 
#' @description
#' Method for calculating production when the PellaTomlinson object is used.
#' 
#' @param biomass FLQuant object.
#' @param params PellaTomlinson object.
#' 
#' @return FLQuant object.
#' 
#' @export
setMethod("production", signature(biomass = "FLQuant", params = "PellaTomlinson"), 
          function(biomass, params) {
            # Use the provided biomass and params from the object
            params_obj = params@params
            
            # Check if we have multiple iterations
            dims = dim(params_obj)
            n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
            
            if (n_iter == 1) {
              # Single iteration - use existing logic
              r = c(params_obj["r"])
              p = c(params_obj["p"])
              virgin = c(params_obj["virgin"])
              return(((r / p) * biomass) * (1 - exp(log(biomass / virgin) * p)))
            } else {
              # Multiple iterations - handle each iteration properly
              # Extract all parameters for all iterations
              r_vals = c(params_obj["r"])
              p_vals = c(params_obj["p"])
              virgin_vals = c(params_obj["virgin"])
              
              result_list = list()
              for (i in 1:n_iter) {
                r_i = r_vals[i]
                p_i = p_vals[i]
                virgin_i = virgin_vals[i]
                prod_i = ((r_i / p_i) * biomass) * (1 - exp(log(biomass / virgin_i) * p_i))
                result_list[[i]] = prod_i
              }
              # Combine results into a single FLQuant with iterations
              proto = result_list[[1]]
              dims_proto = dim(proto)
              dn_proto = dimnames(proto)
              dims_proto[6] = length(result_list)
              dn_proto$iter = as.character(seq_len(length(result_list)))
              result = FLQuant(array(NA, dim = dims_proto, dimnames = dn_proto))
              for (i in seq_along(result_list)) {
                result[, , , , , i] = result_list[[i]]
              }
              return(result)
            }
          })

#' Method for MSY (PellaTomlinson)
#' 
#' @description
#' Method for calculating MSY when the PellaTomlinson object is used.
#' 
#' @param params PellaTomlinson object.
#' 
#' @return Numeric value or FLQuant object.
#' 
#' @export
setMethod("msy", signature(params = "PellaTomlinson"), 
          function(params) {
            # Extract parameters and estimate MSY from r, p, and virgin
            params_obj = params@params
            
            # Check if we have multiple iterations
            dims = dim(params_obj)
            n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
            
            if (n_iter == 1) {
              # Single iteration - return as FLPar for consistency
              r = c(params_obj["r"])
              p = c(params_obj["p"])
              virgin = c(params_obj["virgin"])
              # For Pella-Tomlinson: MSY = (r/p) * BMSY * (1 - (BMSY/virgin)^p)
              bmsy = virgin * exp(log(1 / (1 + p)) / p)
              msy_val = ((r / p) * bmsy) * (1 - (bmsy / virgin)^p)
              return(FLPar(msy = msy_val))
            } else {
              # Multiple iterations - return multi-iteration results
              r_vals = c(params_obj["r"])
              p_vals = c(params_obj["p"])
              virgin_vals = c(params_obj["virgin"])
              
              # Calculate MSY for each iteration
              msy_vals = numeric(n_iter)
              for (i in 1:n_iter) {
                r_i = r_vals[i]
                p_i = p_vals[i]
                virgin_i = virgin_vals[i]
                bmsy_i = virgin_i * exp(log(1 / (1 + p_i)) / p_i)
                msy_vals[i] = ((r_i / p_i) * bmsy_i) * (1 - (bmsy_i / virgin_i)^p_i)
              }
              
              # Return as FLPar with iterations
              result_matrix = matrix(msy_vals, nrow = 1, ncol = n_iter)
              return(FLPar(result_matrix, dimnames = list(params = "msy", iter = 1:n_iter)))
            }
          })

#' Method for BMSY (PellaTomlinson)
#' 
#' @description
#' Method for calculating BMSY when the PellaTomlinson object is used.
#' 
#' @param params PellaTomlinson object.
#' 
#' @return Numeric value or FLQuant object.
#' 
#' @export
setMethod("bmsy", signature(params = "PellaTomlinson"), 
          function(params) {
            # Extract parameters and estimate BMSY from p and virgin
            params_obj = params@params
            
            # Check if we have multiple iterations
            dims = dim(params_obj)
            n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
            
            if (n_iter == 1) {
              # Single iteration - return as FLPar for consistency
              p = c(params_obj["p"])
              virgin = c(params_obj["virgin"])
              # For Pella-Tomlinson: BMSY = virgin * (1/(1+p))^(1/p)
              bmsy_val = virgin * exp(log(1 / (1 + p)) / p)
              return(FLPar(bmsy = bmsy_val))
            } else {
              # Multiple iterations - return multi-iteration results
              p_vals = c(params_obj["p"])
              virgin_vals = c(params_obj["virgin"])
              
              # Calculate BMSY for each iteration
              bmsy_vals = numeric(n_iter)
              for (i in 1:n_iter) {
                p_i = p_vals[i]
                virgin_i = virgin_vals[i]
                bmsy_vals[i] = virgin_i * exp(log(1 / (1 + p_i)) / p_i)
              }
              
              # Return as FLPar with iterations
              result_matrix = matrix(bmsy_vals, nrow = 1, ncol = n_iter)
              return(FLPar(result_matrix, dimnames = list(params = "bmsy", iter = 1:n_iter)))
            }
          })

#' Method for FMSY (PellaTomlinson)
#' 
#' @description
#' Method for calculating FMSY when the PellaTomlinson object is used.
#' 
#' @param params PellaTomlinson object.
#' 
#' @return Numeric value or FLQuant object.
#' 
#' @export
setMethod("fmsy", signature(params = "PellaTomlinson"), 
          function(params) {
            # Extract parameters and estimate FMSY from r and p
            params_obj = params@params
            
            # Check if we have multiple iterations
            dims = dim(params_obj)
            n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
            
            if (n_iter == 1) {
              # Single iteration - return as FLPar for consistency
              r = c(params_obj["r"])
              p = c(params_obj["p"])
              fmsy_val = r / p
              return(FLPar(fmsy = fmsy_val))
            } else {
              # Multiple iterations - return multi-iteration results
              r_vals = c(params_obj["r"])
              p_vals = c(params_obj["p"])
              
              # Calculate FMSY for each iteration
              fmsy_vals = numeric(n_iter)
              for (i in 1:n_iter) {
                r_i = r_vals[i]
                p_i = p_vals[i]
                fmsy_vals[i] = r_i / p_i
              }
              
              # Return as FLPar with iterations
              result_matrix = matrix(fmsy_vals, nrow = 1, ncol = n_iter)
              return(FLPar(result_matrix, dimnames = list(params = "fmsy", iter = 1:n_iter)))
            }
          })

#' p method for PellaTomlinson class - derives p from bmsy and virgin
#'
#' @description
#' Derives the shape parameter p from BMSY and virgin biomass.
#'
#' @param params PellaTomlinson object.
#'
#' @return FLPar object with p parameter.
#'
#' @export
setMethod("p", signature(params = "PellaTomlinson"),
          function(params) {
            params_obj = params@params
            
            # Check if we have multiple iterations
            dims = dim(params_obj)
            n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
            
            if (n_iter == 1) {
              # Single iteration - return as FLPar for consistency
              bmsy = c(params_obj["bmsy"])
              virgin = c(params_obj["virgin"])
              params_numeric = c(bmsy = bmsy, virgin = virgin)
              p_val = pFn(params_numeric)
              return(FLPar(p = p_val))
            } else {
              # Multiple iterations - return multi-iteration results
              bmsy_vals = c(params_obj["bmsy"])
              virgin_vals = c(params_obj["virgin"])
              
              # Calculate p for each iteration using pFn
              p_vals = numeric(n_iter)
              for (i in 1:n_iter) {
                params_i = c(bmsy = bmsy_vals[i], virgin = virgin_vals[i])
                p_vals[i] = pFn(params_i)
              }
              
              # Return as FLPar with iterations
              result_matrix = matrix(p_vals, nrow = 1, ncol = n_iter)
              return(FLPar(result_matrix, dimnames = list(params = "p", iter = 1:n_iter)))
            }
          })

#' r method for PellaTomlinson class - derives r from msy, bmsy, and virgin
#' 
#' @description
#' Derives the intrinsic growth rate r from MSY, BMSY, and virgin biomass.
#' 
#' @param params PellaTomlinson object.
#' 
#' @return FLPar object with r parameter.
#' 
#' @export
setMethod("r", signature(params = "PellaTomlinson"), 
          function(params) {
            params_obj = params@params
            
            # Check if we have multiple iterations
            dims = dim(params_obj)
            n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
            
                        if (n_iter == 1) {
              # Single iteration - return as FLPar for consistency
              msy = params_obj["msy"]
              bmsy = params_obj["bmsy"]
              virgin = params_obj["virgin"]
              
              # Use pFn utility function to avoid circular dependency
              p_val = pFn(c(bmsy = bmsy, virgin = virgin))
              if (is.na(p_val)) return(FLPar(r = NA))
              
              # For the Pella-Tomlinson model, solve for r directly using FLR operators
              r_val = ((msy %*% p_val) %/% bmsy) %/% (1 - exp(log(bmsy %/% virgin) %*% p_val))
              
              return(FLPar(r = r_val))
            } else {
              # Multiple iterations - return multi-iteration results
              msy_vals = c(params_obj["msy"])
              bmsy_vals = c(params_obj["bmsy"])
              virgin_vals = c(params_obj["virgin"])
              
              # Get p values for all iterations using pFn utility function
              p_vals = numeric(n_iter)
              for (i in 1:n_iter) {
                params_i = c(bmsy = bmsy_vals[i], virgin = virgin_vals[i])
                p_vals[i] = pFn(params_i)
              }
              
              # Calculate r for each iteration
              r_vals = numeric(n_iter)
              for (i in 1:n_iter) {
                if (is.na(p_vals[i])) {
                  r_vals[i] = NA
                  next
                }
                
                msy_i = msy_vals[i]
                bmsy_i = bmsy_vals[i]
                virgin_i = virgin_vals[i]
                p_i = p_vals[i]
                
                # For the Pella-Tomlinson model, solve for r directly
                r_vals[i] = ((msy_i * p_i) / bmsy_i) / (1 - (bmsy_i / virgin_i)^p_i)
              }
              
              # Return as FLPar with iterations
              result_matrix = matrix(r_vals, nrow = 1, ncol = n_iter)
              return(FLPar(result_matrix, dimnames = list(params = "r", iter = 1:n_iter)))
            }
          })

#' Utility function for estimating Pella-Tomlinson parameters from bmsy, virgin, and msy
#' 
#' @description
#' Estimates Pella-Tomlinson parameters (r, p) from observed BMSY, virgin, and MSY.
#' This function is designed to work with ddply and other data manipulation functions.
#' 
#' @param bmsy Numeric vector of BMSY values.
#' @param virgin Numeric vector of virgin biomass values.
#' @param msy Numeric vector of MSY values.
#' 
#' @return Data frame with estimated r and p values.
#' 
#' @examples
#' bmsy_vals = c(0.7, 0.5, 0.8, 0.3)
#' virgin_vals = c(1.0, 1.0, 1.0, 1.0)
#' msy_vals = c(0.1, 0.08, 0.12, 0.06)
#' 
#' estimated_params = estimatePellaTomlinsonParams(bmsy_vals, virgin_vals, msy_vals)
#' 
#' @export
estimatePellaTomlinsonParams = function(bmsy, virgin, msy) {
  # Input validation
  if (length(bmsy) != length(virgin) || length(virgin) != length(msy)) {
    stop("bmsy, virgin, and msy must have the same length")
  }
  
  if (any(bmsy <= 0) || any(virgin <= 0) || any(msy <= 0)) {
    stop("bmsy, virgin, and msy must all be positive")
  }
  
  if (any(bmsy >= virgin)) {
    stop("bmsy must be less than virgin for all observations")
  }
  
  # Initialize results vectors
  n = length(bmsy)
  r_vals = numeric(n)
  p_vals = numeric(n)
  
  # Process each observation
  for (i in 1:n) {
    # Create parameter vector for this observation
    params_i = c(bmsy = bmsy[i], virgin = virgin[i], msy = msy[i])
    
    # Estimate p first
    p_i = pFn(params_i)
    p_vals[i] = p_i
    
    # Then estimate r
    if (!is.na(p_i)) {
      r_i = rFn(params_i)
      r_vals[i] = r_i
    } else {
      r_vals[i] = NA
    }
  }
  
  # Return as data frame suitable for ddply
  return(data.frame(
    r = r_vals,
    p = p_vals,
    bmsy = bmsy,
    virgin = virgin,
    msy = msy
  ))
}

#' Alternative version that returns just r and p (for cases where you only need the parameters)
#' 
#' @description
#' Estimates Pella-Tomlinson parameters (r, p) from observed BMSY, virgin, and MSY.
#' This function is designed to work with ddply and other data manipulation functions.
#' 
#' @param bmsy Numeric vector of BMSY values.
#' @param virgin Numeric vector of virgin biomass values.
#' @param msy Numeric vector of MSY values.
#' 
#' @return Data frame with estimated r and p values.
#' 
#' @examples
#' bmsy_vals = c(0.7, 0.5, 0.8, 0.3)
#' virgin_vals = c(1.0, 1.0, 1.0, 1.0)
#' msy_vals = c(0.1, 0.08, 0.12, 0.06)
#' 
#' estimated_params = estimatePellaTomlinsonParamsSimple(bmsy_vals, virgin_vals, msy_vals)
#' 
#' @export
estimatePellaTomlinsonParamsSimple = function(bmsy, virgin, msy) {
  # Input validation
  if (length(bmsy) != length(virgin) || length(virgin) != length(msy)) {
    stop("bmsy, virgin, and msy must have the same length")
  }
  
  if (any(bmsy <= 0) || any(virgin <= 0) || any(msy <= 0)) {
    stop("bmsy, virgin, and msy must all be positive")
  }
  
  if (any(bmsy >= virgin)) {
    stop("bmsy must be less than virgin for all observations")
  }
  
  # Initialize results vectors
  n = length(bmsy)
  r_vals = numeric(n)
  p_vals = numeric(n)
  
  # Process each observation
  for (i in 1:n) {
    # Create parameter vector for this observation
    params_i = c(bmsy = bmsy[i], virgin = virgin[i], msy = msy[i])
    
    # Estimate p first
    p_i = pFn(params_i)
    p_vals[i] = p_i
    
    # Then estimate r
    if (!is.na(p_i)) {
      r_i = rFn(params_i)
      r_vals[i] = r_i
    } else {
      r_vals[i] = NA
    }
  }
  
  # Return just r and p
  return(data.frame(r = r_vals, p = p_vals))
}

#' Fast parameter extraction utility
#' 
#' @description
#' Quickly derives parameters from MSY, BMSY, and virgin biomass and returns them as numeric.
#' Vectorised over inputs. If the recycled length is 1, returns a named numeric vector; otherwise,
#' returns a data.frame with columns r, p, virgin, bmsy, msy.
#' 
#' @param msy Numeric vector or FLPar. Maximum Sustainable Yield.
#' @param bmsy Numeric vector or FLPar. Biomass at MSY.
#' @param virgin Numeric vector or FLPar. Virgin biomass/carrying capacity.
#' 
#' @return If length 1 after recycling: named numeric vector (r, p, virgin, bmsy, msy).
#' Otherwise: data.frame with columns (r, p, virgin, bmsy, msy).
#' 
#' @examples
#' # Single set (returns named numeric vector)
#' fastParams(msy = 0.1, bmsy = 0.3, virgin = 1)
#' 
#' # Vectorised over virgin (returns data.frame)
#' fastParams(msy = 0.1, bmsy = 0.3, virgin = seq(1, 2, 0.1))
#' 
#' # With FLPar object
#' fastParams(FLPar(msy = 0.1, bmsy = 0.3, virgin = seq(1, 2, 0.1)))
#' 
#' @export
fastParams = function(msy, bmsy, virgin) {
  # Determine intended length via recycling rules
  len = max(length(msy), length(bmsy), length(virgin))
  msy_vec = rep(msy, length.out = len)
  bmsy_vec = rep(bmsy, length.out = len)
  virgin_vec = rep(virgin, length.out = len)

  # Basic validation
  if (any(msy_vec <= 0) || any(bmsy_vec <= 0) || any(virgin_vec <= 0)) {
    stop("msy, bmsy, and virgin must all be positive")
  }
  if (any(bmsy_vec >= virgin_vec)) {
    stop("bmsy must be less than virgin")
  }

  # Compute p and r for each set
  p_vals = vapply(seq_len(len), function(i) {
    params_i = c(bmsy = bmsy_vec[i], virgin = virgin_vec[i], msy = msy_vec[i])
    unname(as.numeric(pFn(params_i)))
  }, numeric(1))

  r_vals = vapply(seq_len(len), function(i) {
    params_i = c(bmsy = bmsy_vec[i], virgin = virgin_vec[i], msy = msy_vec[i])
    unname(as.numeric(rFn(params_i)))
  }, numeric(1))

  if (len == 1) {
    out = c(r_vals[1], p_vals[1], virgin_vec[1], bmsy_vec[1], msy_vec[1])
    names(out) = c("r", "p", "virgin", "bmsy", "msy")
    return(out)
  }

  res = data.frame(r = r_vals, p = p_vals, virgin = virgin_vec, bmsy = bmsy_vec, msy = msy_vec)
  return(res)
}

#' @rdname fastParams
#' @export
setMethod("fastParams", signature(msy = "FLPar", bmsy = "missing", virgin = "missing"),
          function(msy, bmsy, virgin) {
            # Extract parameters from FLPar
            param_names = dimnames(msy)$params
            
            if (!all(c("msy", "bmsy", "virgin") %in% param_names)) {
              stop("FLPar must contain: msy, bmsy, virgin")
            }
            
            # Extract values
            msy_vals = c(msy["msy"])
            bmsy_vals = c(msy["bmsy"])
            virgin_vals = c(msy["virgin"])
            
            # Call the numeric version
            fastParams(msy = msy_vals, bmsy = bmsy_vals, virgin = virgin_vals)
          })

#' Utility function to create data frame with stock biomass and yield values
#' 
#' @description
#' Creates a data frame with biomass and yield values for plotting and analysis.
#' 
#' @param params Numeric vector containing 'r', 'p', and 'virgin' parameters.
#' @param n_steps Number of biomass steps (default: 101).
#' 
#' @return Data frame with biomass and yield values.
#' 
#' @examples
#' params = c(r = 0.5, p = 0.25, virgin = 1.0)
#' data_frame = createPellaTomlinsonData(params)
#' 
#' @export
createPellaTomlinsonData = function(..., n_steps = 101) {
  # Input validation
  args = list(...)
  
  # Handle both FLPar and individual parameters
  if (length(args) == 1 && inherits(args[[1]], "FLPar")) {
    # Single FLPar argument
    params = args[[1]]
    required_params = c("r", "p", "virgin")
    if (!all(required_params %in% dimnames(params)$params)) {
      stop("params must contain: r, p, virgin")
    }
    r = c(params["r"])
    p = c(params["p"])
    virgin = c(params["virgin"])
  } else if (length(args) >= 3) {
    # Individual parameters
    if (length(args) == 3) {
      r = args[[1]]
      p = args[[2]]
      virgin = args[[3]]
    } else {
      # Try to extract by name if more than 3 arguments
      if ("r" %in% names(args) && "p" %in% names(args) && "virgin" %in% names(args)) {
        r = args[["r"]]
        p = args[["p"]]
        virgin = args[["virgin"]]
      } else {
        stop("Must provide either an FLPar object or r, p, and virgin parameters")
      }
    }
  } else {
    stop("Must provide either an FLPar object or r, p, and virgin parameters")
  }
  
  if (is.na(r) || is.na(p) || is.na(virgin)) {
    stop("r, p, and virgin must not be NA")
  }
  
  if (virgin <= 0) {
    stop("virgin must be positive")
  }
  
  if (n_steps < 2) {
    stop("n_steps must be at least 2")
  }
  
  # Create biomass sequence from 0 to virgin
  biomass = seq(0, virgin, length.out = n_steps)
  
  # Calculate yield/production for each biomass level
  # Create a simple FLPar for the production calculation
  params_flpar = FLPar(r = r, p = p, virgin = virgin)
  biomass_flquant = FLQuant(biomass)
  yield = production(biomass_flquant, params_flpar)
  
  # Create data frame
  result = data.frame(
    biomass = biomass,
    yield = c(yield)
  )
  
  # Add additional useful columns
  result$biomass_ratio = result$biomass / virgin
  result$yield_per_biomass = ifelse(result$biomass > 0, result$yield / result$biomass, 0)
  
  # Add MSY and BMSY reference points
  msy_val = msy(params_flpar)
  bmsy_val = bmsy(params_flpar)
  
  # Add flags for reference points
  result$is_bmsy = abs(result$biomass - bmsy_val) < (virgin / (n_steps - 1)) / 2
  result$is_msy = abs(result$yield - msy_val) < (max(yield, na.rm = TRUE) / (n_steps - 1)) / 2
  
  return(result)
}

#' Print method for PellaTomlinson class
#' 
#' @description
#' Prints a summary of the PellaTomlinson object.
#' 
#' @param x PellaTomlinson object.
#' @param ... Additional arguments.
#' 
#' @export
setMethod("print", signature(x = "PellaTomlinson"),
          function(x, ...) {
            cat("Pella-Tomlinson Production Model\n")
            cat("================================\n")
            cat("Parameters:\n")
            print(x@params)
            cat("\nStock biomass range:", range(x@stock), "\n")
            cat("Stock biomass dimensions:", dim(x@stock), "\n")
          })

#' Show method for PellaTomlinson class
#' 
#' @description
#' Shows the summary of the PellaTomlinson object.
#' 
#' @param object PellaTomlinson object.
#' 
#' @export
setMethod("show", signature(object = "PellaTomlinson"),
          function(object) {
            print(object)
          })

#' Summary method for PellaTomlinson class
#' 
#' @description
#' Provides a summary of the PellaTomlinson object, including model parameters,
#' derived reference points, and stock data.
#' 
#' @param object PellaTomlinson object.
#' @param ... Additional arguments.
#' 
#' @export
setMethod("summary", signature(object = "PellaTomlinson"),
          function(object, ...) {
            cat("Pella-Tomlinson Production Model Summary\n")
            cat("========================================\n")
            cat("Model Parameters:\n")
            cat("  r (intrinsic growth rate):", round(c(object@params["r"]), 6), "\n")
            cat("  p (shape parameter):", round(c(object@params["p"]), 6), "\n")
            cat("  virgin (carrying capacity):", round(c(object@params["virgin"]), 6), "\n")
            cat("  bmsy (biomass at MSY):", round(c(object@params["bmsy"]), 6), "\n")
            cat("  msy (maximum sustainable yield):", round(c(object@params["msy"]), 6), "\n")
            cat("\nDerived Reference Points:\n")
            cat("  MSY (calculated):", round(c(msy(object)), 6), "\n")
            cat("  BMSY (calculated):", round(c(bmsy(object)), 6), "\n")
            cat("  FMSY (calculated):", round(c(fmsy(object)), 6), "\n")
            cat("\nStock Data:\n")
            cat("  Biomass range:", round(range(object@stock), 4), "\n")
            cat("  Number of points:", length(object@stock), "\n")
            cat("  Dimensions:", paste(dim(object@stock), collapse = " x "), "\n")
          })

#' Plot method for PellaTomlinson class using ggplot
#' 
#' @description
#' Plots the production curve of the PellaTomlinson object.
#' 
#' @param x PellaTomlinson object.
#' @param y Missing (production is plotted against biomass).
#' @param ... Additional arguments for ggplot.
#' 
#' @export
setMethod("plot", signature(x = "PellaTomlinson", y = "missing"),
          function(x, y, ...) {
            # Check if ggplot2 is available
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
              stop("ggplot2 package is required for plotting. Please install it.")
            }
            
            # Check if we have multiple iterations
            dims = dim(x@params)
            n_iter = ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
            
            if (n_iter == 1) {
              # Single iteration - use existing logic
              # Generate production curve
              biomass = c(x@stock)
              production_vals = c(production(x))
              
              # Create data frame for plotting
              plot_data = data.frame(
                biomass = biomass,
                production = production_vals
              )
              
              # Calculate reference points for annotation
              msy_point = c(msy(x))
              bmsy_point = c(bmsy(x))
              fmsy_point = c(fmsy(x))
              
              # Create the plot
              p = ggplot2::ggplot(plot_data, ggplot2::aes(x = biomass, y = production)) +
                ggplot2::geom_line(color = "steelblue", size = 1) +
                ggplot2::geom_point(data = data.frame(x = bmsy_point, y = msy_point), 
                                   ggplot2::aes(x = x, y = y), 
                                   color = "red", size = 3, shape = 17) +
                ggplot2::geom_vline(xintercept = bmsy_point, linetype = "dashed", color = "red", alpha = 0.7) +
                ggplot2::geom_hline(yintercept = msy_point, linetype = "dashed", color = "red", alpha = 0.7) +
                ggplot2::annotate("text", x = bmsy_point * 1.1, y = msy_point * 1.05,
                                 label = paste("MSY =", round(msy_point, 4)), 
                                 color = "red", size = 4) +
                ggplot2::annotate("text", x = bmsy_point * 1.1, y = msy_point * 0.95,
                                 label = paste("BMSY =", round(bmsy_point, 4)), 
                                 color = "red", size = 4) +
                ggplot2::labs(
                  title = "Pella-Tomlinson Production Curve",
                  subtitle = paste("r =", round(c(x@params["r"]), 4), 
                                 ", p =", round(c(x@params["p"]), 4),
                                 ", virgin =", round(c(x@params["virgin"]), 4)),
                  x = "Biomass",
                  y = "Production",
                  caption = paste("FMSY =", round(fmsy_point, 4))
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                  plot.title = ggplot2::element_text(size = 14, face = "bold"),
                  plot.subtitle = ggplot2::element_text(size = 10),
                  axis.title = ggplot2::element_text(size = 12),
                  axis.text = ggplot2::element_text(size = 10)
                )
              
              return(p)
            } else {
              # Multiple iterations - plot first iteration for now
              # This could be enhanced to show multiple curves
              warning("Multiple iterations detected. Plotting first iteration only.")
              
              # Use first iteration parameters
              biomass = c(x@stock)
              production_vals = c(production(x))
              
              # Create data frame for plotting
              plot_data = data.frame(
                biomass = biomass,
                production = production_vals
              )
              
              # Calculate reference points for annotation (first iteration)
              msy_point = c(msy(x))
              bmsy_point = c(bmsy(x))
              fmsy_point = c(fmsy(x))
              
              # Create the plot
              p = ggplot2::ggplot(plot_data, ggplot2::aes(x = biomass, y = production)) +
                ggplot2::geom_line(color = "steelblue", size = 1) +
                ggplot2::geom_point(data = data.frame(x = bmsy_point, y = msy_point), 
                                   ggplot2::aes(x = x, y = y), 
                                   color = "red", size = 3, shape = 17) +
                ggplot2::geom_vline(xintercept = bmsy_point, linetype = "dashed", color = "red", alpha = 0.7) +
                ggplot2::geom_hline(yintercept = msy_point, linetype = "dashed", color = "red", alpha = 0.7) +
                ggplot2::annotate("text", x = bmsy_point * 1.1, y = msy_point * 1.05,
                                 label = paste("MSY =", round(msy_point, 4)), 
                                 color = "red", size = 4) +
                ggplot2::annotate("text", x = bmsy_point * 1.1, y = msy_point * 0.95,
                                 label = paste("BMSY =", round(bmsy_point, 4)), 
                                 color = "red", size = 4) +
                ggplot2::labs(
                  title = "Pella-Tomlinson Production Curve (First Iteration)",
                  subtitle = paste("r =", round(c(x@params["r"]), 4), 
                                 ", p =", round(c(x@params["p"]), 4),
                                 ", virgin =", round(c(x@params["virgin"]), 4),
                                 "\nMultiple iterations available: ", n_iter),
                  x = "Biomass",
                  y = "Production",
                  caption = paste("FMSY =", round(fmsy_point, 4))
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                  plot.title = ggplot2::element_text(size = 14, face = "bold"),
                  plot.subtitle = ggplot2::element_text(size = 10),
                  axis.title = ggplot2::element_text(size = 12),
                  axis.text = ggplot2::element_text(size = 10)
                )
              
              return(p)
            }
          })

#' Forward projection method for PellaTomlinson class
#' 
#' @description
#' Projects stock biomass year-by-year for given fishing mortality (F) or catch.
#' 
#' @param object PellaTomlinson object.
#' @param F_val Fishing mortality rate (optional, if not provided, catch must be).
#' @param catch Catch level (optional, if not provided, F must be).
#' @param years Number of years to project (default: 50).
#' @param initial_biomass Initial biomass (default: virgin biomass).
#' @param ... Additional arguments.
#' 
#' @return FLQuant object with biomass projections over time.
#' 
#' @details
#' The forward projection uses the Pella-Tomlinson production function:
#' B(t+1) = B(t) + P(B(t)) - C(t)
#' where P(B) is production and C(t) is catch.
#' 
#' @examples
#' pt_model = PellaTomlinson(FLPar(msy = 0.1, bmsy = 0.4, virgin = 1))
#' 
#' # Project with constant F
#' biomass_proj = fwd(pt_model, F_val = 0.2, years = 30)
#' 
#' # Project with constant catch
#' biomass_proj = fwd(pt_model, catch = 0.05, years = 30)


#' @rdname fwd
#' @export
setMethod("fwd", signature(object = "PellaTomlinson"),
          function(object, F_val = NULL, catch = NULL, years = 50, initial_biomass = NULL, ...) {
            # Extract parameters
            params_obj = object@params
            r_val = c(params_obj["r"])
            p_val = c(params_obj["p"])
            virgin = c(params_obj["virgin"])
            
            # Set initial biomass
            if (is.null(initial_biomass)) {
              initial_biomass = virgin
            }
            
            # Validate inputs
            if (is.null(F_val) && is.null(catch)) {
              stop("Either F_val or catch must be provided")
            }
            
            if (years < 1) {
              stop("years must be at least 1")
            }
            
            if (initial_biomass <= 0) {
              stop("initial_biomass must be positive")
            }
            
            # Initialize biomass vector
            biomass = numeric(years + 1)
            biomass[1] = initial_biomass
            
            # Project biomass year by year
            for (t in 1:years) {
              # Current biomass
              B_t = biomass[t]
              
              # Calculate production
              if (abs(p_val) < 1e-10) {
                # Fox model approximation for p near 0
                production_t = r_val * B_t * (1 - B_t / virgin)
              } else {
                production_t = (r_val / p_val) * B_t * (1 - (B_t / virgin)^p_val)
              }
              
              # Calculate catch
              if (!is.null(F_val)) {
                # F-based catch
                catch_t = F_val * B_t
              } else {
                # Constant catch
                catch_t = catch
              }
              
              # Ensure catch doesn't exceed available biomass + production
              max_catch = max(0, B_t + production_t)
              catch_t = min(catch_t, max_catch)
              
              # Update biomass for next year
              biomass[t + 1] = max(0, B_t + production_t - catch_t)
            }
            
            # Create FLQuant with proper dimensions
            result = FLQuant(biomass, 
                           dimnames = list(year = 0:years, 
                                         iter = 1))
            
            return(result)
          })

setGeneric("recoveryPT", function(object, ...) standardGeneric("recoveryPT"))

#' @rdname recoveryPT
#' @export
setMethod("recoveryPT", signature(object = "PellaTomlinson"),
          function(object, F_val, initial_biomass = NULL, 
                   target_biomass = NULL, max_years = 100, 
                   tolerance = 0.01, ...) {
            # Extract parameters
            params_obj = object@params
            r_val = c(params_obj["r"])
            p_val = c(params_obj["p"])
            virgin = c(params_obj["virgin"])
            bmsy_val = c(bmsy(object))
            
            # Set default values
            if (is.null(initial_biomass)) {
              initial_biomass = virgin * 0.1  # Start at 10% of virgin biomass
            }
            
            if (is.null(target_biomass)) {
              target_biomass = bmsy_val
            }
            
            # Validate inputs
            if (F_val < 0) {
              stop("F_val must be non-negative")
            }
            
            if (initial_biomass <= 0) {
              stop("initial_biomass must be positive")
            }
            
            if (target_biomass <= 0) {
              stop("target_biomass must be positive")
            }
            
            if (max_years < 1) {
              stop("max_years must be at least 1")
            }
            
            if (tolerance <= 0) {
              stop("tolerance must be positive")
            }
            
            # Initialize variables
            biomass = numeric(max_years + 1)
            biomass[1] = initial_biomass
            years_to_recovery = NA
            reached_target = FALSE
            
            # Project forward until target is reached or max_years exceeded
            for (t in 1:max_years) {
              # Current biomass
              B_t = biomass[t]
              
              # Calculate production
              if (abs(p_val) < 1e-10) {
                # Fox model approximation for p near 0
                production_t = r_val * B_t * (1 - B_t / virgin)
              } else {
                production_t = (r_val / p_val) * B_t * (1 - (B_t / virgin)^p_val)
              }
              
              # Calculate catch
              catch_t = F_val * B_t
              
              # Ensure catch doesn't exceed available biomass + production
              max_catch = max(0, B_t + production_t)
              catch_t = min(catch_t, max_catch)
              
              # Update biomass for next year
              biomass[t + 1] = max(0, B_t + production_t - catch_t)
              
              # Check if target is reached
              if (abs(biomass[t + 1] - target_biomass) <= tolerance) {
                years_to_recovery = t
                reached_target = TRUE
                break
              }
              
              # Check if biomass is stable (no change)
              if (abs(biomass[t + 1] - B_t) < 1e-6) {
                break
              }
            }
            
            # Create FLQuant with biomass trajectory
            biomass_trajectory = FLQuant(biomass[1:(t + 1)], 
                                         dimnames = list(year = 0:t, 
                                                       iter = 1))
            
            # Return results
            result = list(
              years = years_to_recovery,
              biomass_trajectory = biomass_trajectory,
              reached_target = reached_target,
              final_biomass = biomass[t + 1]
            )
            
            return(result)
          })

#' @rdname recoveryPT
#' @export
setMethod("recoveryPT", signature(object = "numeric"),
          function(object, ...) {
            # Get additional arguments passed via ...
            additional_args <- list(...)
            
            # For numeric vectors, ensure parameter names are correct
            if ("b0" %in% names(object)) {
              # Convert b0 to initial for compatibility with recovery function
              names(object)[names(object) == "b0"] <- "initial"
            }
            
            # Combine object with additional arguments
            if (length(additional_args) > 0) {
              # Create a combined parameter list to preserve vector nature
              combined_params <- as.list(object)
              # Add additional arguments
              for (name in names(additional_args)) {
                if (name %in% names(combined_params)) {
                  # If parameter already exists, replace it
                  combined_params[[name]] <- additional_args[[name]]
                } else {
                  # If parameter doesn't exist, add it
                  combined_params[[name]] <- additional_args[[name]]
                }
              }
            } else {
              combined_params <- as.list(object)
            }
            
            # Check if we can use the optimized version
            # Look for scalar bmsy, msy, virgin but vector initial and/or ftar
            scalar_params <- c("bmsy", "msy", "virgin", "r")
            vector_params <- c("initial", "ftar")
            
            # Add default values for missing vector parameters
            if (!("initial" %in% names(combined_params))) {
              combined_params$initial <- 0.1
            }
            if (!("ftar" %in% names(combined_params))) {
              combined_params$ftar <- 0
            }
            
            param_names <- names(combined_params)
            has_scalar_vector_mix <- FALSE
            
            # Check if we have scalar parameters
            has_scalars <- any(sapply(scalar_params, function(p) p %in% param_names))
            
            # Check if we have vector parameters
            has_vectors <- any(sapply(vector_params, function(p) p %in% param_names && length(combined_params[[p]]) > 1))
            
            # If we have both scalar and vector parameters, use optimized version
            if (has_scalars && has_vectors) {
              # Check that scalar parameters are actually single values
              scalar_vals <- sapply(scalar_params, function(p) if(p %in% param_names) length(combined_params[[p]]) == 1 else TRUE)
              if (all(scalar_vals)) {
                has_scalar_vector_mix <- TRUE
              }
            }
        
            if (has_scalar_vector_mix) {
              print("Use optimized version for better performance")
              return(recovery_optimized(combined_params))
            } else {
              print("Use standard recovery function")
              return(recovery(combined_params))
            }
          })

#' @rdname recoveryPT
#' @export
setMethod("recoveryPT", signature(object = "data.frame"),
          function(object, ...) {
            # Get additional arguments passed via ...
            additional_args <- list(...)
            
            # For data.frames, ensure parameter names are correct
            if ("b0" %in% names(object)) {
              # Convert b0 to initial for compatibility with recovery function
              names(object)[names(object) == "b0"] <- "initial"
            }
            
            # Combine object with additional arguments
            if (length(additional_args) > 0) {
              # Create a combined parameter list
              combined_params <- as.list(object)
              # Add additional arguments
              for (name in names(additional_args)) {
                if (name %in% names(combined_params)) {
                  # If parameter already exists, replace it
                  combined_params[[name]] <- additional_args[[name]]
                } else {
                  # If parameter doesn't exist, add it
                  combined_params[[name]] <- additional_args[[name]]
                }
              }
            } else {
              combined_params <- as.list(object)
            }
            
            # Check if we can use the optimized version
            # Look for columns that should be scalar vs vector
            scalar_params <- c("bmsy", "msy", "virgin", "r")
            vector_params <- c("initial", "ftar")
            
            param_names <- names(combined_params)
            has_scalar_vector_mix <- FALSE
            
            # Check if we have scalar parameters (should be single unique values)
            has_scalars <- any(sapply(scalar_params, function(p) p %in% param_names))
            
            # Check if we have vector parameters (should have multiple unique values)
            has_vectors <- any(sapply(vector_params, function(p) p %in% param_names && length(combined_params[[p]]) > 1))
            
            # If we have both scalar and vector parameters, use optimized version
            if (has_scalars && has_vectors) {
              # Check that scalar parameters are actually single values
              scalar_vals <- sapply(scalar_params, function(p) {
                if(p %in% param_names) {
                  # For list, check if it's a single value
                  val <- combined_params[[p]]
                  return(length(val) == 1)
                } else {
                  return(TRUE)
                }
              })
              
              if (all(scalar_vals)) {
                has_scalar_vector_mix <- TRUE
              }
            }
            
            if (has_scalar_vector_mix) {
              # Use optimized version for better performance
              return(recovery_optimized(combined_params))
            } else {
              # Use standard recovery function
              return(recovery(combined_params))
            }
          })

#' @rdname recoveryPT
#' @export
setMethod("recoveryPT", signature(object = "FLPar"),
          function(object, ...) {
            # Get additional arguments passed via ...
            additional_args <- list(...)
            
            # For FLPar objects, ensure parameter names are correct
            param_names <- dimnames(object)$params
            if ("b0" %in% param_names) {
              # Convert b0 to initial for compatibility with recovery function
              param_names[param_names == "b0"] <- "initial"
              dimnames(object)$params <- param_names
            }
            
            # If there are additional arguments, we need to handle them
            # For FLPar, we'll convert to a list and combine with additional args
            if (length(additional_args) > 0) {
              # Convert FLPar to list with proper names
              param_names <- dimnames(object)$params
              params_list <- list()
              for (i in seq_along(param_names)) {
                params_list[[param_names[i]]] <- as.numeric(object[i])
              }
              # Add additional arguments
              for (name in names(additional_args)) {
                if (name %in% names(params_list)) {
                  # If parameter already exists, replace it
                  params_list[[name]] <- additional_args[[name]]
                } else {
                  # If parameter doesn't exist, add it
                  params_list[[name]] <- additional_args[[name]]
                }
              }
              
              # Check if we can use the optimized version
              # Look for scalar bmsy, msy, virgin but vector initial and/or ftar
              scalar_params <- c("bmsy", "msy", "virgin", "r")
              vector_params <- c("initial", "ftar")
              
              param_names <- names(params_list)
              has_scalar_vector_mix <- FALSE
              
              # Debug: Print what we're checking
              cat("Debug FLPar method:\n")
              cat("  param_names:", paste(param_names, collapse = ", "), "\n")
              cat("  scalar_params:", paste(scalar_params, collapse = ", "), "\n")
              cat("  vector_params:", paste(vector_params, collapse = ", "), "\n")
              
              # Check if we have scalar parameters
              has_scalars <- any(sapply(scalar_params, function(p) p %in% param_names))
              cat("  has_scalars:", has_scalars, "\n")
              
              # Check if we have vector parameters
              has_vectors <- any(sapply(vector_params, function(p) p %in% param_names && length(params_list[[p]]) > 1))
              cat("  has_vectors:", has_vectors, "\n")
              
              # If we have both scalar and vector parameters, use optimized version
              if (has_scalars && has_vectors) {
                # Check that scalar parameters are actually single values
                scalar_vals <- sapply(scalar_params, function(p) if(p %in% param_names) length(params_list[[p]]) == 1 else TRUE)
                cat("  scalar_vals:", paste(scalar_vals, collapse = ", "), "\n")
                if (all(scalar_vals)) {
                  has_scalar_vector_mix <- TRUE
                }
              }
              
              cat("  has_scalar_vector_mix:", has_scalar_vector_mix, "\n")
              
              if (has_scalar_vector_mix) {
                # Use optimized version for better performance
                cat("  Using recovery_optimized\n")
                return(recovery_optimized(params_list))
              } else {
                # Use standard recovery function
                cat("  Using recovery\n")
                return(recovery(params_list))
              }
            } else {
              # No additional arguments, use FLPar directly
              return(recovery(object))
            }
          })

#' @rdname recoveryPT
#' @export
setMethod("recoveryPT", signature(object = "list"),
          function(object, ...) {
            # Get additional arguments passed via ...
            additional_args <- list(...)
            
            # For lists, ensure parameter names are correct
            if ("b0" %in% names(object)) {
              # Convert b0 to initial for compatibility with recovery function
              names(object)[names(object) == "b0"] <- "initial"
            }
            
            # Combine object with additional arguments
            if (length(additional_args) > 0) {
              # Create a combined parameter list
              combined_params <- object
              # Add additional arguments
              for (name in names(additional_args)) {
                if (name %in% names(combined_params)) {
                  # If parameter already exists, replace it
                  combined_params[[name]] <- additional_args[[name]]
                } else {
                  # If parameter doesn't exist, add it
                  combined_params[[name]] <- additional_args[[name]]
                }
              }
            } else {
              combined_params <- object
            }
            
            # Check if we can use the optimized version
            # Look for scalar bmsy, msy, virgin but vector initial and/or ftar
            scalar_params <- c("bmsy", "msy", "virgin", "r")
            vector_params <- c("initial", "ftar")
            
            param_names <- names(combined_params)
            has_scalar_vector_mix <- FALSE
            
            # Check if we have scalar parameters
            has_scalars <- any(sapply(scalar_params, function(p) p %in% param_names))
            
            # Check if we have vector parameters
            has_vectors <- any(sapply(vector_params, function(p) p %in% param_names && length(combined_params[[p]]) > 1))
            
            # If we have both scalar and vector parameters, use optimized version
            if (has_scalars && has_vectors) {
              # Check that scalar parameters are actually single values
              scalar_vals <- sapply(scalar_params, function(p) {
                if(p %in% param_names) {
                  return(length(combined_params[[p]]) == 1)
                } else {
                  return(TRUE)
                }
              })
              
              if (all(scalar_vals)) {
                has_scalar_vector_mix <- TRUE
              }
            }
            
            if (has_scalar_vector_mix) {
              # Use optimized version for better performance
              return(recovery_optimized(combined_params))
            } else {
              # Use standard recovery function
              return(recovery(combined_params))
            }
          })






# Recovery function for calculating time to reach target biomass
# This function calculates how long it takes for a stock to recover to BMSY under given fishing mortality

.recovery_single <- function(
    params,
    virgin = 1
) {
  # Extract parameters from the named vector
  r <- params["r"]
  initial <- params["initial"]
  ftar <- params["ftar"]
  
  # If virgin is in params, use it; otherwise use the default
  if ("virgin" %in% names(params)) {
    virgin <- params["virgin"]
  }
  
  # Handle msy input - convert to r if needed
  if (is.na(r) && "msy" %in% names(params)) {
    msy_val <- params["msy"]
    if (is.finite(msy_val)) {
      # We need bmsy and p to calculate r from msy
      # First check if we have bmsy or p
      val_bmsy <- if ("bmsy" %in% names(params)) as.numeric(params["bmsy"]) else NA_real_
      val_p    <- if ("p"    %in% names(params)) as.numeric(params["p"])    else NA_real_
      
      if (is.finite(val_bmsy)) {
        # We have bmsy, estimate p first, then r
        bmsy <- val_bmsy
        target_ratio <- bmsy / virgin
        if (!is.finite(target_ratio) || target_ratio <= 0 || target_ratio >= 1) {
          stop("bmsy must be a fraction in (0,1) of virgin (K)")
        }
        
        fn <- function(p) (1/(1+p))^(1/p) - target_ratio
        p_est <- tryCatch({
          uniroot(fn, interval=c(1e-6, 100), tol=1e-8, extendInt="upX")$root
        }, error=function(e) {
          warning("Numerical solution for p failed, using default p = 0.25")
          0.25
        })
        
        # Now calculate r from msy: r = msy * p / (bmsy * (1 - (bmsy/virgin)^p))
        r <- msy_val * p_est / (bmsy * (1 - (bmsy / virgin)^p_est))
      } else if (is.finite(val_p)) {
        # We have p, calculate bmsy, then r
        p_est <- val_p
        bmsy <- virgin * (1/(1+p_est))^(1/p_est)
        r <- msy_val * p_est / (bmsy * (1 - (bmsy / virgin)^p_est))
      } else {
        stop("Must provide either 'bmsy' or 'p' along with 'msy' to calculate 'r'")
      }
    }
  }
  
  # Check which parameters are provided and estimate accordingly
  val_bmsy <- if ("bmsy" %in% names(params)) as.numeric(params["bmsy"]) else NA_real_
  val_p    <- if ("p"    %in% names(params)) as.numeric(params["p"])    else NA_real_
  has_bmsy <- is.finite(val_bmsy)
  has_p    <- is.finite(val_p)
  
  if (has_bmsy && has_p) {
    # Both provided - use them directly
    bmsy <- val_bmsy
    p_est <- val_p
  } else if (has_bmsy && !has_p) {
    # bmsy provided, estimate p from BMSY/K = (1/(1+p))^(1/p)
    bmsy <- params["bmsy"]
    target_ratio <- bmsy / virgin
    if (!is.finite(target_ratio) || target_ratio <= 0 || target_ratio >= 1) {
      stop("bmsy must be a fraction in (0,1) of virgin (K)")
    }
    
    fn <- function(p) (1/(1+p))^(1/p) - target_ratio
    p_est <- tryCatch({
      uniroot(fn, interval=c(1e-6, 100), tol=1e-8, extendInt="upX")$root
    }, error=function(e) {
      warning("Numerical solution for p failed, using default p = 0.25")
      0.25
    })
  } else if (has_p && !has_bmsy) {
    # p provided, estimate bmsy from BMSY = virgin * (1/(1+p))^(1/p)
    p_est <- val_p
    bmsy <- virgin * (1/(1+p_est))^(1/p_est)
  } else {
    # Neither provided - error
    stop("Must provide either 'bmsy' or 'p' in params")
  }

  # Check that r is valid
  if (is.na(r) || !is.finite(r)) {
    stop("Could not calculate 'r' from the provided parameters. Need either 'r' directly or 'msy' with 'bmsy' or 'p'.")
  }
  
  # Derived reference points
  Bmsy <- virgin * (1/(1+p_est))^(1/p_est)
  Fmsy <- r / (1 + p_est)  # Correct: Fmsy = MSY/Bmsy = r/(1+p)
  F <- ftar * Fmsy
  B0 <- initial * virgin
  # Compute MSY consistently from derived points
  msy_val <- Fmsy * Bmsy

  # --- Step 2: set up Pella-Tomlinson ODE ---
  # dB/dt = (r/p) * B * (1 - (B/K)^p) - F * B
  if (!requireNamespace("deSolve", quietly=TRUE)) stop("Package 'deSolve' is required")

  pellatomlinson <- function(t, state, parms) {
    B <- state[1]
    dB <- (r/p_est) * B * (1 - (B/virgin)^p_est) - F * B
    list(c(dB))
  }

  times <- seq(0, 500, by = 0.1)
  out <- deSolve::ode(y=c(B=B0), times=times, func=pellatomlinson, parms=NULL)

  # --- Step 3: find time to reach BMSY ---
  B_traj <- out[, 2]  # Second column is the biomass values
  hit_idx <- which(B_traj >= Bmsy)
  T_recover <- if (length(hit_idx) > 0) out[hit_idx[1], 1] else NA_real_  # First column is time

  return(setNames(
    c(T_recover, msy_val, Fmsy, bmsy, virgin, r, p_est, F, initial),
    c("T_recover", "msy", "Fmsy", "bmsy", "virgin", "r", "p", "F", "initial")
  ))
}

recovery <- function(
    params,   # named numeric vector OR data.frame/matrix with cols: r, b0, ftar, (bmsy|p), optional virgin
    virgin = 1
) {
  # FLPar input (returns FLPar)
  if (inherits(params, "FLPar")) {
    if (!requireNamespace("FLCore", quietly = TRUE)) {
      stop("Package 'FLCore' is required to use FLPar inputs")
    }
    # extract iterations
    dims <- dim(params)
    n_iter <- ifelse(length(dims) > 1 && dims[2] > 1, dims[2], 1)
         # extract vectors per parameter (may be NULL if not present)
     get_vec <- function(nm) if (nm %in% dimnames(params)$params) c(params[nm]) else rep(NA_real_, n_iter)
     r_vec     <- get_vec("r")
     msy_vec   <- get_vec("msy")
     initial_vec    <- get_vec("initial")
     ftar_vec  <- get_vec("ftar")
    bmsy_vec  <- get_vec("bmsy")
    p_vec     <- get_vec("p")
    virgin_vec<- if ("virgin" %in% dimnames(params)$params) c(params["virgin"]) else rep(virgin, n_iter)

         results <- vector("list", n_iter)
     for (i in seq_len(n_iter)) {
       prm_i <- c(r = r_vec[i], initial = initial_vec[i], ftar = ftar_vec[i])
       if (is.finite(msy_vec[i]))   prm_i <- c(prm_i, msy = msy_vec[i])
      if (is.finite(bmsy_vec[i])) prm_i <- c(prm_i, bmsy = bmsy_vec[i])
      if (is.finite(p_vec[i]))    prm_i <- c(prm_i, p = p_vec[i])
      if (is.finite(virgin_vec[i])) prm_i <- c(prm_i, virgin = virgin_vec[i])
      results[[i]] <- .recovery_single(prm_i, virgin = virgin_vec[i])
    }
    # Always return FLPar when input is FLPar (consistent with FLPar input)
    # Multiple iterations: build FLPar with params along first dim and iterations along second
    mat <- do.call(cbind, results)  # rows are named metrics, cols are iters
    if (is.null(dim(mat))) {
      mat <- matrix(mat, ncol = 1)
      rownames(mat) <- names(results[[1]])
      colnames(mat) <- "iter1"
    } else {
      # Ensure matrix has proper column names
      if (is.null(colnames(mat))) {
        colnames(mat) <- paste0("iter", 1:ncol(mat))
      }
    }
    
    # Ensure matrix has correct structure for FLPar
    # FLPar expects params as first dimension, iterations as second dimension
    # The issue is that FLPar constructor adds an extra dimension when given a matrix
    # We need to ensure the matrix is properly structured
    if (ncol(mat) == 1) {
      # Single iteration - create FLPar with correct dimensions
      # For single iteration, we need to ensure the matrix is 2D, not 1D
      # The problem is that FLPar constructor is interpreting the matrix dimensions incorrectly
      # We need to create it with explicit dimensions
      flp <- FLCore::FLPar(array(mat, dim = c(nrow(mat), 1)),
        dimnames = list(params = rownames(mat), iter = "1"))
    } else {
      # Multiple iterations
      flp <- FLCore::FLPar(array(mat, dim = c(nrow(mat), ncol(mat))),
        dimnames = list(params = rownames(mat), iter = seq_len(ncol(mat))))
    }
    return(flp)
  }

  # Vectorized entry points:
  # - data.frame/matrix: rows are scenarios
  # - list of vectors: build data.frame by recycling to max length
  # - named numeric vector: single scenario

  # data.frame or matrix of scenarios
  if (is.data.frame(params) || (is.matrix(params) && !is.null(colnames(params)))) {
    prm_df <- as.data.frame(params, stringsAsFactors = FALSE)
    nrows <- nrow(prm_df)
    results <- vector("list", nrows)
    for (i in seq_len(nrows)) {
      row_vec <- as.numeric(prm_df[i, , drop = TRUE])
      names(row_vec) <- colnames(prm_df)
      vrg <- if ("virgin" %in% names(row_vec)) row_vec["virgin"] else virgin
      results[[i]] <- .recovery_single(row_vec, virgin = vrg)
    }
    out <- do.call(rbind, lapply(results, function(x) as.data.frame(as.list(x))))
    rownames(out) <- rownames(prm_df)
    return(out)
  }

  # list of vectors (possibly different lengths) -> build data.frame by recycling
  if (is.list(params) && !is.null(names(params)) && is.null(dim(params))) {
    lens <- vapply(params, length, integer(1))
    maxlen <- if (length(lens) > 0) max(lens) else 0
    if (maxlen > 1) {
      prm_df <- as.data.frame(lapply(params, function(x) rep(x, length.out = maxlen)), stringsAsFactors = FALSE)
      return(recovery(prm_df, virgin = virgin))
    }
    # fall-through to single scenario
  }

  # single named numeric vector
  if (is.numeric(params) && is.null(dim(params))) {
    return(.recovery_single(params, virgin = if ("virgin" %in% names(params)) params["virgin"] else virgin))
  }

  stop("Unsupported 'params' type. Use a named numeric vector, a data.frame/matrix with columns, or a named list of vectors.")
}

# Example usage (numeric vector):
# params = c(r=0.1, bmsy=0.4, b0=0.1, ftar=0.1)
# result = recovery(params)
# print(result)
# 
# # Or with p instead of bmsy:
# params = c(r=0.1, p=0.2, b0=0.1, ftar=0.1)
# result = recovery(params)
# print(result)
# 
# # Or with explicit virgin:
# params = c(r=0.1, bmsy=0.4, b0=0.1, ftar=0.1, virgin=2)
# result = recovery(params)
# print(result)

#' Optimized Recovery Function for Multiple initial/ftar Values
#' 
#' @description
#' Optimized version of recovery function that estimates p and r only once when
#' bmsy, msy, and virgin are single values, but initial and/or ftar are vectors.
#' This is much more efficient than calling recovery() multiple times.
#' 
#' @param params Named list or data.frame where:
#'   - bmsy, msy, virgin: single numeric values (scalars)
#'   - initial, ftar: numeric vectors of any length
#'   - r: optional single value (if not provided, will be estimated from msy, bmsy, virgin)
#' @param virgin Default virgin biomass if not in params (default: 1)
#' 
#' @return Data frame with recovery results for all combinations of initial and ftar
#' 
#' @examples
#' # Single bmsy, msy, virgin but multiple initial and ftar values
#' params <- list(
#'   bmsy = 0.4,      # single value
#'   msy = 0.1,       # single value  
#'   virgin = 2,      # single value
#'   initial = c(0.1, 0.2, 0.3),  # vector
#'   ftar = c(0.1, 0.2)           # vector
#' )
#' result <- recovery_optimized(params)
#' 
#' # With data.frame input
#' df <- data.frame(
#'   bmsy = 0.4, msy = 0.1, virgin = 2,  # these will be recycled
#'   initial = c(0.1, 0.2, 0.3),
#'   ftar = c(0.1, 0.2)
#' )
#' result <- recovery_optimized(df)
#' 
#' @export
recovery_optimized <- function(params, virgin = 1) {
  # Convert to list if data.frame
  if (is.data.frame(params)) {
    params <- as.list(params)
  }
  
  # Extract scalar parameters (should be single values)
  scalar_params <- c("bmsy", "msy", "virgin", "r")
  vector_params <- c("initial", "ftar")
  
  # Check which parameters are provided
  param_names <- names(params)
  
  # Validate scalar parameters are single values
  for (param in scalar_params) {
    if (param %in% param_names) {
      if (length(params[[param]]) != 1) {
        stop(paste("Parameter", param, "must be a single value, not a vector"))
      }
    }
  }
  
  # Check that we have at least one vector parameter
  has_vector <- any(sapply(vector_params, function(p) p %in% param_names && length(params[[p]]) > 1))
  if (!has_vector) {
    warning("No vector parameters found. Consider using recovery() instead.")
    return(recovery(params, virgin = virgin))
  }
  
  # Extract scalar values
  bmsy_val <- if ("bmsy" %in% param_names) params$bmsy else NA_real_
  msy_val <- if ("msy" %in% param_names) params$msy else NA_real_
  virgin_val <- if ("virgin" %in% param_names) params$virgin else virgin
  r_val <- if ("r" %in% param_names) params$r else NA_real_
  
  # Extract vector parameters
  initial_vals <- if ("initial" %in% param_names) params$initial else 0.1
  ftar_vals <- if ("ftar" %in% param_names) params$ftar else 0
  
  # Create combinations only if both are vectors, otherwise use single values
  if (length(initial_vals) == 1 && length(ftar_vals) == 1) {
    # Both are single values - single scenario
    combinations <- data.frame(
      initial = initial_vals,
      ftar = ftar_vals,
      stringsAsFactors = FALSE
    )
  } else if (length(initial_vals) == 1) {
    # Only ftar is a vector
    combinations <- data.frame(
      initial = rep(initial_vals, length(ftar_vals)),
      ftar = ftar_vals,
      stringsAsFactors = FALSE
    )
  } else if (length(ftar_vals) == 1) {
    # Only initial is a vector
    combinations <- data.frame(
      initial = initial_vals,
      ftar = rep(ftar_vals, length(initial_vals)),
      stringsAsFactors = FALSE
    )
  } else {
    # Both are vectors - create all combinations
    combinations <- expand.grid(
      initial = initial_vals,
      ftar = ftar_vals,
      stringsAsFactors = FALSE
    )
  }
  
  n_combinations <- nrow(combinations)
  
  # Estimate p and r only once if needed
  if (is.na(r_val)) {
    if (is.na(bmsy_val) || is.na(msy_val)) {
      stop("Must provide either 'r' directly or both 'bmsy' and 'msy' to estimate r")
    }
    
    # Estimate p from bmsy and virgin
    target_ratio <- bmsy_val / virgin_val
    if (target_ratio <= 0 || target_ratio >= 1) {
      stop("bmsy must be a fraction in (0,1) of virgin")
    }
    
    fn <- function(p) (1/(1+p))^(1/p) - target_ratio
    p_est <- tryCatch({
      uniroot(fn, interval=c(1e-6, 100), tol=1e-8, extendInt="upX")$root
    }, error=function(e) {
      warning("Numerical solution for p failed, using default p = 0.25")
      0.25
    })
    # Calculate r from msy, bmsy, virgin, and p
    r_val <- msy_val * p_est / (bmsy_val * (1 - (bmsy_val / virgin_val)^p_est))
  } else {
    # r is provided, estimate p from bmsy and virgin
    if (is.na(bmsy_val)) {
      stop("Must provide 'bmsy' to estimate p when r is provided")
    }
    
    target_ratio <- bmsy_val / virgin_val
    if (target_ratio <= 0 || target_ratio >= 1) {
      stop("bmsy must be a fraction in (0,1) of virgin")
    }
    
    fn <- function(p) (1/(1+p))^(1/p) - target_ratio
    p_est <- tryCatch({
      uniroot(fn, interval=c(1e-6, 100), tol=1e-8, extendInt="upX")$root
    }, error=function(e) {
      warning("Numerical solution for p failed, using default p = 0.25")
      0.25
    })
  }
  
  # Derived reference points (calculated once)
  Bmsy <- virgin_val * (1/(1+p_est))^(1/p_est)
  Fmsy <- r_val / (1 + p_est)
  msy_calculated <- Fmsy * Bmsy
  
  # Calculate recovery for each combination
  results <- vector("list", n_combinations)
  
  for (i in 1:n_combinations) {
    initial_i <- combinations$initial[i]
    ftar_i <- combinations$ftar[i]
    
    # Calculate F for this combination
    F_i <- ftar_i * Fmsy
    # initial should be a fraction of BMSY, not virgin biomass
    B0_i <- initial_i * Bmsy
    
    # Set up ODE for this combination
    if (!requireNamespace("deSolve", quietly=TRUE)) {
      stop("Package 'deSolve' is required")
    }
    
    pellatomlinson <- function(t, state, parms) {
      B <- state[1]
      dB <- (r_val/p_est) * B * (1 - (B/virgin_val)^p_est) - F_i * B
      list(c(dB))
    }
    
    times <- seq(0, 500, by = 0.1)
    out <- deSolve::ode(y=c(B=B0_i), times=times, func=pellatomlinson, parms=NULL)
    
    # Find time to reach BMSY
    B_traj <- out[, 2]
    hit_idx <- which(B_traj >= Bmsy)
    T_recover <- if (length(hit_idx) > 0) out[hit_idx[1], 1] else NA_real_
    
    if (is.na(T_recover)) {
      warning("No recovery time found for initial=", initial_i, ", ftar=", ftar_i, 
              ". B0=", B0_i, ", BMSY=", Bmsy, ", max B=", max(B_traj))
    }
    
    # Store results
    results[[i]] <- list(
      T_recover = T_recover,
      msy = msy_calculated,
      Fmsy = Fmsy,
      bmsy = Bmsy,
      virgin = virgin_val,
      r = r_val,
      p = p_est,
      F = F_i,
      initial = initial_i
    )
    

  }
  
  # Combine results into data frame
  result_df <- do.call(rbind, lapply(results, function(x) as.data.frame(as.list(x))))
  
  # Add ftar column for clarity
  result_df$ftar <- combinations$ftar
  
  # Ensure all expected columns exist
  expected_cols <- c("T_recover", "msy", "Fmsy", "bmsy", "virgin", "r", "p", "F", "initial", "ftar")
  missing_cols <- setdiff(expected_cols, names(result_df))
  
  if (length(missing_cols) > 0) {
    warning("Missing columns in result: ", paste(missing_cols, collapse = ", "))
  }
  
  # Reorder columns for consistency with recovery(), keeping only existing columns
  existing_cols <- names(result_df)
  col_order <- expected_cols[expected_cols %in% existing_cols]
  result_df <- result_df[, col_order]
  
  return(result_df)
}

#' Extract Reference Points from FLBRP Object
#' 
#' @description Extracts key reference points (MSY, crash, virgin) from an FLBRP object
#' and returns them as an FLPar object with harvest, yield, and SSB values.
#' 
#' @param x An FLBRP object
#' @return FLPar object containing reference points
#' @export
#' @examples
#' \dontrun{
#' # Assuming you have an FLBRP object called 'brp'
#' refs(brp)
#' }
setGeneric("refs", function(x,...) standardGeneric("refs"))

#' @rdname refs
#' @export
setMethod("refs", signature(x="FLBRP"), function(x,shape=NULL) {
  rfs=refpts(x)[c("msy", "crash", "virgin"), c("harvest", "yield", "ssb")]
  
  rtn=rbind(msy   =rfs[1, "yield"],
            bmsy  =rfs[1, "ssb"],
            virgin=rfs[3, "ssb"],
            fcrash=rfs[2, "harvest"],
            fmsy  =rfs[1, "harvest"])
  
  if (!is.null(shape)) rtn["virgin"]=rtn["bmsy"]/shape
  
  FLPar(rtn[drop=TRUE])})

setGeneric("refsEB", function(x,...) standardGeneric("refsEB"))

#' @rdname refsEB
#' @export
setMethod("refsEB", signature(x="FLBRP"), function(x,shape=NULL) {
  rfs=refptsEB(x)[c("msy", "crash", "virgin"), c("harvest", "yield", "eb")]
  
  rtn=rbind(msy   =rfs[1, "yield"],
            bmsy  =rfs[1, "eb"],
            virgin=rfs[3, "eb"],
            fcrash=rfs[2, "harvest"],
            fmsy  =rfs[1, "harvest"])
  
  if (!is.null(shape)) rtn["virgin"]=rtn["bmsy"]/shape
  
  FLPar(rtn[drop=TRUE])})

