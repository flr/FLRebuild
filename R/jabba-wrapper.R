auxFn <- function(lag = 0, obsE = 0.3, sigma = TRUE, type = "", ...) {
  args <- list(...)
   
  auxiliary <- NULL
  
  if (any(names(args) %in% c("z", "f", "ffmsy", "effort", "bbmsy", "bk"))) {
    type <- names(args)[names(args) %in% c("z", "f", "ffmsy", "effort", "bbmsy", "bk")][1]
    auxiliary <- args[[type]]
  }
  
  list(
    auxiliary = auxiliary,
    auxiliary.type = type,
    auxiliary.sigma = sigma,  # estimated?
    auxiliary.obsE = obsE,
    auxiliary.lag = lag  # lag effect between impact and Z pop structure
  )
}

mFn <- function(shape, fmsy) {
  mi <- seq(0.01, 2, 0.001)
  m <- (mi^(-1/(mi-1)) - shape)^2
  m <- mi[m == min(m)]
  r <- (1 - exp(-fmsy)) * (m - 1) / (1 - m^-1)
  c(m = m, r = r)
}

#' Calculate JABBA Prior Parameters from ICES Data
#'
#' @description
#' Extracts and calculates prior parameters for JABBA (Just Another Bayesian
#' Biomass Assessment) models from ICES stock assessment data. This function
#' combines information from equilibrium simulations, benchmark calculations,
#' and FishLife predictions to create a comprehensive set of priors.
#'
#' @param icesdata A list of ICES stock assessment objects (typically FLStock
#'   or similar objects) with names corresponding to stock IDs
#'
#' @return A data.frame with prior parameters for each stock, containing:
#'   \itemize{
#'     \item \code{.id}: Stock identifier
#'     \item \code{r}: Intrinsic growth rate (from FishLife)
#'     \item \code{shape}: BMSY/B0 ratio (shape parameter for production function)
#'     \item \code{psi}: Initial depletion (initial SSB / B0)
#'     \item \code{ssb.maxyear}: Current depletion (current SSB / BMSY)
#'     \item \code{ssb.minyear}: Initial depletion relative to BMSY
#'     \item \code{initial}: Initial SSB value
#'     \item \code{current}: Current SSB value
#'     \item Additional columns from benchmark calculations (e.g., fmsy, bmsy, b0)
#'   }
#'
#' @details
#' This function integrates multiple data sources:
#' \itemize{
#'   \item \strong{Equilibrium simulations:} Provides BMSY and B0 values
#'   \item \strong{Benchmark calculations:} Provides reference points (FMSY, etc.)
#'   \item \strong{FishLife:} Provides prior estimates for intrinsic growth rate (r)
#'   \item \strong{Time series:} Extracts initial and current SSB values
#' }
#'
#' The shape parameter (BMSY/B0) is a key input for the Pella-Tomlinson production
#' function used in JABBA models.
#'
#' @examples
#' \dontrun{
#' # Calculate priors for multiple stocks
#' priors <- jabbaPriors(icesdata)
#' 
#' # Use priors in JABBA model fitting
#' fit <- jabba(catch = catch_data, pr = priors[priors$.id == "cod.27.47d20", ])
#' }
#'
#' @seealso \code{\link{jabba}}, \code{\link{jabbaData}}
#'
#' @export
jabbaPriors <- function(icesdata) {
  eqsm <- eqsim(icesdata)
  benchm <- benchmark(icesdata)
  initial <- ldply(icesdata, function(x) {
    rtn <- tseries(x)
    rtn[rtn$year == min(rtn$year), ]
  })
  current <- ldply(icesdata, function(x) {
    rtn <- tseries(x)
    rtn[rtn$year == max(rtn$year), ]
  })
  fl <- fishlife(icesdata)
  
  priors <- merge(benchm[, -8], eqsm[, c(".id", "bmsy", "b0")], by = ".id")
  priors <- merge(priors, transmute(fl, .id = .id, r = r), by = ".id")
  priors <- merge(priors, transmute(initial, .id = .id, initial = ssb), by = ".id")
  priors <- merge(priors, transmute(current, .id = .id, current = ssb), by = ".id")
  priors <- transform(priors,
                      shape = bmsy / b0,
                      ssb.maxyear = current / bmsy,
                      ssb.minyear = initial / b0,
                      psi = initial / b0)
  
  priors
}

#' Prepare Data for JABBA Model Fitting
#'
#' @description
#' Extracts and formats time series data from ICES stock assessment objects
#' for use in JABBA (Just Another Bayesian Biomass Assessment) models.
#' This function prepares catch, biomass index, and fishing mortality data,
#' with optional support for historical catch data and survey indices.
#'
#' @param id Character string. Stock identifier matching a name in \code{icesdata}
#' @param icesdata A list of ICES stock assessment objects (typically FLStock
#'   or similar objects) with names corresponding to stock IDs
#' @param ctc1903 Optional data.frame. Historical catch data (1903 onwards)
#'   with columns: \code{.id}, \code{year}, and \code{catch}. If provided,
#'   replaces standard catch data for years where available.
#' @param indices Optional data.frame. Survey indices with columns: \code{.id},
#'   \code{year}, \code{survey}, and \code{data}. Multiple surveys are supported
#'   and will be cast to wide format (one column per survey).
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{year}: Year
#'     \item \code{catch}: Catch data (missing or zero values replaced with small positive value)
#'     \item \code{index}: Estimated biomass index (ebiomass)
#'     \item \code{ffmsy}: Fishing mortality relative to FMSY
#'     \item Additional columns for each survey (if \code{indices} provided)
#'   }
#'
#' @details
#' The function performs several data preparation steps:
#' \itemize{
#'   \item Extracts catch time series and replaces missing/zero values with
#'     a small positive value (1e-6 times mean catch)
#'   \item Extracts estimated biomass index (ebiomass) from time series
#'   \item Calculates F/FMSY from fishing mortality and FMSY reference point
#'   \item Optionally merges historical catch data (1903+) if provided
#'   \item Optionally adds survey indices, casting to wide format with one
#'     column per survey name
#' }
#'
#' Missing catch values are replaced to ensure JABBA models can run, as JABBA
#' requires positive catch values for all years in the time series.
#'
#' @examples
#' \dontrun{
#' # Basic usage with ICES data
#' data <- jabbaData("cod.27.47d20", icesdata)
#' 
#' # With historical catch data
#' data <- jabbaData("cod.27.47d20", icesdata, ctc1903 = historical_catch)
#' 
#' # With survey indices
#' data <- jabbaData("cod.27.47d20", icesdata, indices = survey_data)
#' 
#' # Full example
#' data <- jabbaData("cod.27.47d20", icesdata, 
#'                   ctc1903 = historical_catch, 
#'                   indices = survey_data)
#' }
#'
#' @seealso \code{\link{jabba}}, \code{\link{jabbaPriors}}
#'
#' @export
jabbaData <- function(id, icesdata, ctc1903 = NULL, indices = NULL) {
  ts <- model.frame(tseries(icesdata[[id]]), drop = TRUE)
  
  catch <- ts[, c("year", "catch")]
  names(catch)[2] <- "catch"
  small <- mean(ts$catch, na.rm = TRUE) * 1e-6
  catch$catch[is.na(catch$catch) | catch$catch <= 0] <- small
  
  eb <- transmute(ts, year = year, index = eb)
  eb[is.na(eb)] <- NA
  
  fmsy <- benchmark(icesdata[[id]])["fmsy"]
  ffmsy <- transmute(ts, year = year, ffmsy = (1 - exp(-f)) / (1 - exp(-c(fmsy))))
  
  rtn <- merge(merge(catch, eb, by = "year"), ffmsy, by = "year")
  
  if (!is.null(ctc1903)) {
    ctc1903 <- subset(ctc1903, .id == id & year <= dims(icesdata[[id]])$maxyear)[, 2:3]
    if (dim(ctc1903)[1] > 0) {
      rtn <- merge(rtn[, -2], ctc1903, by = "year", all.y = TRUE)
      rtn$catch[is.na(rtn$catch) | rtn$catch <= 0] <- small
    }
  }
  
  if (is.null(indices)) return(rtn)
  
  idx <- subset(indices, .id == id)
  
  if (dim(idx)[1] > 0) {
    idx <- try(cast(subset(idx, year <= dims(icesdata[[id]])$maxyear),
                    year ~ survey, value = "data", fun = "mean"))
    idx <- merge(idx, catch, by = "year", all.y = TRUE)[, seq(dim(idx)[2])]
    idx[is.na(idx)] <- NA
    
    return(merge(rtn, idx, by = "year", all.x = TRUE))
  }
  
  return(rtn)
}

#' Fit JABBA (Just Another Bayesian Biomass Assessment) Model
#'
#' @description
#' Fits a JABBA (Just Another Bayesian Biomass Assessment) model to catch and
#' biomass index data. This is a wrapper function that prepares JABBA inputs
#' and calls the underlying JABBA fitting functions (\code{build_jabba} and
#' \code{fit_jabba}).
#'
#' @param catch Data.frame with columns \code{year} and \code{catch}. Catch time
#'   series data. Missing or zero values should be replaced with small positive
#'   values (see \code{\link{jabbaData}}).
#' @param pr Named vector or data.frame row. Prior parameters containing at
#'   minimum: \code{r} (intrinsic growth rate), \code{psi} (initial depletion),
#'   and \code{shape} (BMSY/B0). May also include \code{k} (carrying capacity)
#'   and \code{current} (current depletion). Typically from \code{\link{jabbaPriors}}.
#' @param pr.sd Named vector or data.frame row. Standard deviations for prior
#'   parameters. Must have same names as \code{pr}. Default: 30\% CV for all
#'   parameters (0.3 * pr / pr, which equals 0.3).
#' @param model Character string. JABBA model type. Default: \code{"Pella_m"}
#'   (Pella-Tomlinson with estimated shape parameter m). Other options include
#'   \code{"Schaefer"} (m=2), \code{"Fox"} (m=1), etc.
#' @param assessment Character string. Assessment name/identifier for output.
#'   Default: \code{""}
#' @param scenario Character string. Scenario name/identifier for output.
#'   Default: \code{""}
#' @param index Data.frame or matrix. Biomass index (CPUE) data with years as
#'   rows and surveys as columns. First column should be \code{year}. If NULL,
#'   model uses catch-only data. Default: NULL
#' @param q_bounds Optional list. Bounds for catchability (q) parameters for
#'   each survey index. Default: NULL (uses \code{q_bound} for all surveys)
#' @param sigma.est Logical. Whether to estimate observation error (sigma.obs).
#'   If FALSE, uses \code{fixed.obsE}. Default: TRUE
#' @param fixed.obsE Numeric. Fixed observation error standard deviation if
#'   \code{sigma.est = FALSE}. Default: 0.1
#' @param sigma.proc Logical. Whether to estimate process error (sigma.proc).
#'   If FALSE, uses \code{fixed.procE}. Default: TRUE
#' @param fixed.procE Numeric. Fixed process error standard deviation if
#'   \code{sigma.proc = FALSE}. Default: 0.3
#' @param igamma Numeric vector of length 2. Inverse-gamma prior parameters
#'   for process error variance: c(shape, rate). Default: c(3.0, 0.1)
#' @param Found Numeric vector of length 2. Bounds for catchability (q)
#'   parameters: c(lower, upper). Default: c(1e-3, 1e+3)
#' @param currentDepletion Character string. Type of current depletion prior:
#'   \code{"bbmsy"} (B/BMSY) or \code{"ffmsy"} (F/FMSY). If empty string,
#'   no current depletion prior is used. Default: \code{""}
#' @param initialDepletion Numeric or NA. Initial depletion value (B0/K).
#'   If NA, uses auxiliary data or defaults. Default: NA
#' @param quick Logical. If TRUE, uses quick MCMC settings (fewer iterations).
#'   If FALSE, uses full MCMC. Default: TRUE
#' @param ... Additional arguments passed to \code{auxFn()} for auxiliary data
#'   (e.g., \code{z}, \code{f}, \code{ffmsy}, \code{effort}, \code{bbmsy}, \code{bk})
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{input}: JABBA input object (from \code{build_jabba})
#'     \item \code{fit}: JABBA fit object (from \code{fit_jabba}), containing:
#'       \itemize{
#'         \item \code{pars_posterior}: Posterior samples of model parameters
#'         \item \code{refpts_posterior}: Posterior samples of reference points
#'         \item \code{timeseries}: Time series of biomass, harvest, etc.
#'         \item \code{kobe}: Kobe plot data (B/BMSY vs F/FMSY)
#'         \item \code{kbtrj}: Kobe trajectory
#'       }
#'   }
#'   Returns \code{NULL} if \code{build_jabba} fails, or list with only
#'   \code{input} if \code{fit_jabba} fails.
#'
#' @details
#' This function is a wrapper around the JABBA package functions. It:
#' \itemize{
#'   \item Prepares prior parameters from \code{pr} and \code{pr.sd}
#'   \item Sets up model configuration based on \code{model} type
#'   \item Handles optional auxiliary data (current/initial depletion)
#'   \item Calls \code{build_jabba()} to create JABBA input object
#'   \item Calls \code{fit_jabba()} to fit the model via MCMC
#' }
#'
#' The Pella-Tomlinson production function is used by default, which allows
#' the shape parameter (m) to be estimated. This provides flexibility in
#' modeling different stock production dynamics.
#'
#' @examples
#' \dontrun{
#' # Prepare data and priors
#' data <- jabbaData("cod.27.47d20", icesdata)
#' priors <- jabbaPriors(icesdata)
#' pr <- priors[priors$.id == "cod.27.47d20", ]
#' 
#' # Fit JABBA model
#' fit <- jabba(catch = data[, c("year", "catch")],
#'              pr = pr,
#'              index = data[, c("year", "index")],
#'              scenario = "EBiomass",
#'              assessment = "cod.27.47d20")
#' 
#' # Extract results
#' results <- jabbaExtract(fit)
#' }
#'
#' @seealso
#' \code{\link{jabbaPriors}}, \code{\link{jabbaData}}, \code{\link{jabbaExtract}},
#' \code{\link{jabbaDiagnostics}}, \code{\link{jabbaPE}}
#'
#' @export
jabba <- function(catch,
                  pr,
                  pr.sd = pr / pr * 0.3,
                  model = "Pella_m",
                  assessment = "",
                  scenario = "",
                  index = NULL,
                  q_bounds = NULL,
                  sigma.est = TRUE,
                  fixed.obsE = 0.1,
                  sigma.proc = TRUE,
                  fixed.procE = 0.3,
                  igamma = c(3.0, 0.1),
                  # sigma.est  =FALSE, # for HIST with ICES data
                  # fixed.obsE =0.15, 
                  # sigma.proc =TRUE,
                  # fixed.procE=0.3,
                  # igamma     =c(0.001,0.001),
                  q_bound = c(1e-3, 1e+3),
                  currentDepletion = "",
                  initialDepletion = NA,
                  quick            =TRUE,
                  nc               =3,
                  ...) {
  
  ## priors
  r <- unlist(c(pr[c("r")]))
  r.prior <- c(r, pr.sd["r"])
  
  psi <- unlist(c(pr[c("psi")]))
  if (is.na(psi)) psi <- 0.9
  psi.prior <- c(psi, pr.sd["psi"])
  
  shape <- unlist(c(pr[c("shape")]))
  shape.cv <- pr.sd["shape"]
  
  k.prior <- NA
  if ("k" %in% dimnames(pr)$params) {
    k <- unlist(c(pr["k"]))
    k.prior <- c(k, pr.sd["k"])
  }
  
  if (!is.null(q_bounds)) {
    args <- list(q_bounds = q_bounds)
  } else {
    args <- list()
  }
  
  args <- c(args, list(
    scenario = scenario,
    assessment = assessment,
    model.type = model,
    BmsyK = shape,
    shape.CV = shape.cv,
    catch = catch,
    cpue = index,
    r.prior = r.prior,
    K.prior = k.prior,
    psi.prior = psi.prior,
    sigma.proc = sigma.proc,
    sigma.est = sigma.est,
    fixed.obsE = fixed.obsE,
    igamma = igamma,
    verbose = FALSE
  ))
  args <- args[!is.na(args)]
  
  if (substr(currentDepletion[1], 1, 1) == "b") {
    args <- c(args, list(b.prior = c(c(pr["current"]), pr.sd["current"],
                                      max(catch$year), "bbmsy")))
  }
  if (substr(currentDepletion[1], 1, 1) == "f") {
    args <- c(args, list(b.prior = c(c(pr["current"]), pr.sd["current"],
                                      max(catch$year), "ffmsy")))
  }
  
  if (!is.na(initialDepletion)) {
    # Future implementation for initialDepletion
    aux <- NULL
  } else {
    aux <- auxFn(...)
  }
  
  args <- c(args, aux)


  ## Fit with Catch + Index: Simple Fox with r = Fmsy
  input <- try(do.call("build_jabba", args))
  
  if ("try-error" %in% class(input)) return(NULL)
  
  fit <- try(fit_jabba(input, quickmcmc = quick, verbose = FALSE, nc=nc))
  
  if ("try-error" %in% class(fit)) return(list(input = input))
  
  list(input = input, fit = fit)
}

#' Extract Data from Single JABBA Fit Object
#'
#' @description
#' Internal helper function to extract posterior samples, trajectory, and
#' prior information from a single JABBA fit object. This function handles
#' the extraction logic for individual JABBA fits.
#'
#' @param x JABBA fit object. Can be either:
#'   \itemize{
#'     \item A list with \code{fit} element containing the JABBA model
#'     \item A JABBA fit object directly (with \code{pars_posterior}, etc.)
#'   }
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{posteriors}: data.frame with posterior samples of parameters
#'       and reference points (merged from \code{pars_posterior},
#'       \code{refpts_posterior}, and \code{kobe})
#'     \item \code{trajectory}: data.frame with Kobe trajectory (year, stock,
#'       harvest, etc.)
#'     \item \code{priors}: named vector with prior parameter values (r, K, psi, m)
#'   }
#'   Returns \code{NULL} if input is NULL or extraction fails.
#'
#' @details
#' This function safely extracts data from JABBA fit objects using tryCatch
#' to handle missing components gracefully. It automatically handles both
#' wrapped (list with $fit) and direct JABBA fit object structures.
#'
#' @keywords internal
jabbaExtractFn <- function(x) {
  # Return NULL if x is NULL
  if (is.null(x)) return(NULL)
  
  if ("fit" %in% names(x)) x <- x$fit
  
  # Safely extract posteriors with NULL checking
  posteriors <- tryCatch({
    if (!is.null(x)) {
      cbind(x$pars_posterior,
            x$refpts_posterior,
            x$kobe)
    } else NULL
  }, error = function(e) NULL)
  
  # Safely extract trajectory
  trajectory <- tryCatch({
    if (!is.null(x)) x$kbtrj else NULL
  }, error = function(e) NULL)
  
  if (!is.null(trajectory)) {
    names(trajectory)[names(trajectory) == "yr"] <- "year"
  }
  
  # Safely extract priors
  priors <- tryCatch({
    if (!is.null(x$settings)) {
      prior_vals <- c(unlist(x$settings[c("r.pr", "K.pr", "psi.pr")]),
                      unlist(x$settings[c("mu.m", "m.CV")]))
      names(prior_vals) <- c("r", "r.pr", "k", "k.pr", "psi", "psi.pr", "m", "m.pr")
      prior_vals
    } else NULL
  }, error = function(e) NULL)
  
  list(posteriors = posteriors,
       trajectory = trajectory,
       priors = priors)
}

#' Extract Data from List of JABBA Fit Objects
#'
#' @description
#' Internal helper function to extract and combine data from a list of JABBA
#' fit objects. Each element in the list should be a single JABBA fit object
#' (or list with $fit element). Results are combined into data.frames with
#' an identifier column (.id).
#'
#' @param jabbaList List of JABBA fit objects. Can be named or unnamed.
#'   If unnamed, identifiers are created as character sequence numbers.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{posteriors}: data.frame with posterior samples, including
#'       \code{.id} column for identification
#'     \item \code{trajectory}: data.frame with Kobe trajectory, including
#'       \code{.id} column
#'     \item \code{priors}: data.frame with prior parameters, including
#'       \code{.id} column
#'   }
#'   Components are NULL if no valid data is found.
#'
#' @details
#' This function processes multiple JABBA fits in a single list structure.
#' It uses \code{jabbaExtract()} for each element and combines results using
#' \code{rbindlist()} from the data.table package for efficiency.
#'
#' @keywords internal
jabbaExtractList <- function(jabbaList) {
  library(data.table)
  
  # Handle both named and unnamed lists efficiently
  listNames <- names(jabbaList)
  if (is.null(listNames)) {
    listNames <- as.character(seq_along(jabbaList))
  }
  
  # Single pass: extract and collect components
  posteriorsList <- list()
  trajectoryList <- list()
  priorsList <- list()
  
  for (i in seq_along(jabbaList)) {
    nm <- listNames[i]
    res <- jabbaExtract(jabbaList[[i]])
    if (is.null(res)) next
    
    # Add identifier and collect components in single pass
    if (!is.null(res$posteriors)) {
      posteriorsList[[length(posteriorsList) + 1]] <- 
        cbind(.id = nm, as.data.frame(res$posteriors))
    }
    if (!is.null(res$trajectory)) {
      trajectoryList[[length(trajectoryList) + 1]] <- 
        cbind(.id = nm, as.data.frame(res$trajectory))
    }
    if (!is.null(res$priors)) {
      priorsList[[length(priorsList) + 1]] <- 
        cbind(.id = nm, as.data.frame(t(res$priors)))
    }
  }
  
  # Combine into data frames efficiently
  combined <- list(
    posteriors = if (length(posteriorsList) > 0) {
      rbindlist(posteriorsList, fill = TRUE)
    } else NULL,
    trajectory = if (length(trajectoryList) > 0) {
      rbindlist(trajectoryList, fill = TRUE)
    } else NULL,
    priors = if (length(priorsList) > 0) {
      rbindlist(priorsList, fill = TRUE)
    } else NULL
  )
  
  if (!is.null(combined$trajectory)) {
    names(combined$trajectory)[names(combined$trajectory) == "yr"] <- "year"
  }
  
  return(combined)
}

#' Extract Data from Nested List of JABBA Fit Objects
#'
#' @description
#' Internal helper function to extract and combine data from a nested list
#' structure of JABBA fit objects. The structure is: list(stock = list(scenario = fit)).
#' Results are combined with both stock (.id) and scenario identifiers.
#'
#' @param listOflists Nested list structure where:
#'   \itemize{
#'     \item Outer level: stock identifiers (names)
#'     \item Inner level: scenario identifiers (names) -> JABBA fit objects
#'   }
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{posteriors}: data.frame with posterior samples, including
#'       \code{.id} (stock) and \code{Scenario} columns
#'     \item \code{trajectory}: data.frame with Kobe trajectory, including
#'       \code{.id} and \code{Scenario} columns
#'     \item \code{priors}: data.frame with prior parameters, including
#'       \code{.id} and \code{Scenario} columns
#'   }
#'   Components are NULL if no valid data is found.
#'
#' @details
#' This function handles the most complex nested structure: multiple stocks
#' each with multiple scenarios. It loops through both levels and combines
#' results with appropriate identifiers for downstream analysis.
#'
#' @keywords internal
jabbaExtractLists <- function(listOflists) {
  library(data.table)
  
  # Initialize empty lists for each component
  allPosteriors <- list()
  allTrajectory <- list()
  allPriors <- list()
  
  # Loop through the outer list
  for (.id in names(listOflists)) {
    # Get the inner list
    innerList <- listOflists[[.id]]
    if (is.null(innerList)) next
    
    # Loop through inner lists
    for (Scenario in names(innerList)) {
      # Extract data using your original function
      res <- jabbaExtract(innerList[[Scenario]])
      if (is.null(res)) next
      
      # Add identifiers to non-NULL components
      if (!is.null(res$posteriors)) {
        res$posteriors <- cbind(
          .id = .id,
          Scenario = Scenario,
          as.data.frame(res$posteriors))
        allPosteriors[[paste(.id, Scenario, sep = "_")]] <- res$posteriors
      }
      
      if (!is.null(res$trajectory)) {
        res$trajectory <- cbind(
          .id = .id,
          Scenario = Scenario,
          as.data.frame(res$trajectory))
        allTrajectory[[paste(.id, Scenario, sep = "_")]] <- res$trajectory
      }
      
      if (!is.null(res$priors)) {
        res$priors <- cbind(
          .id = .id,
          Scenario = Scenario,
          as.data.frame(t(res$priors)))
        allPriors[[paste(.id, Scenario, sep = "_")]] <- res$priors
      }
    }
  }
  
  # Combine into final data frames, handling empty lists
  combined <- list(
    posteriors = if (length(allPosteriors) > 0) {
      rbindlist(allPosteriors, fill = TRUE)
    } else NULL,
    trajectory = if (length(allTrajectory) > 0) {
      rbindlist(allTrajectory, fill = TRUE)
    } else NULL,
    priors = if (length(allPriors) > 0) {
      rbindlist(allPriors, fill = TRUE)
    } else NULL
  )
  
  if (!is.null(combined$trajectory)) {
    names(combined$trajectory)[names(combined$trajectory) == "yr"] <- "year"
  }
  
  return(combined)
}

#' Extract Data from JABBA Fit Objects
#'
#' @description
#' Generic function to extract posterior samples, trajectories, and prior
#' information from JABBA (Just Another Bayesian Biomass Assessment) fit
#' objects. This function automatically detects the structure of the input
#' (single fit, list of fits, or nested list) and dispatches to the
#' appropriate extraction method.
#'
#' @param object JABBA fit object(s). Can be:
#'   \itemize{
#'     \item Single JABBA fit object (list with \code{input} and \code{fit},
#'       or direct fit object with \code{pars_posterior})
#'     \item List of JABBA fit objects (for multiple scenarios of one stock)
#'     \item Nested list: list(stock = list(scenario = fit)) (for multiple
#'       stocks each with multiple scenarios)
#'   }
#' @param ... Additional arguments (currently unused)
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{posteriors}: data.frame with posterior samples of parameters
#'       and reference points
#'     \item \code{trajectory}: data.frame with Kobe trajectory (B/BMSY vs F/FMSY)
#'     \item \code{priors}: data.frame or vector with prior parameter values
#'   }
#'   For nested structures, includes identifier columns (\code{.id}, \code{Scenario}).
#'   Returns \code{NULL} if input is NULL or structure cannot be determined.
#'
#' @details
#' This generic function automatically detects the structure of JABBA fit objects:
#' \itemize{
#'   \item \strong{Single fit:} Extracts directly using \code{jabbaExtractFn()}
#'   \item \strong{List of fits:} Extracts each element and combines with
#'     \code{jabbaExtractList()}
#'   \item \strong{Nested list:} Extracts from nested structure with
#'     \code{jabbaExtractLists()}
#' }
#'
#' The function checks for standard JABBA object components:
#' \code{pars_posterior}, \code{refpts_posterior}, \code{kobe}, \code{kbtrj}.
#'
#' @examples
#' \dontrun{
#' # Single JABBA fit
#' fit <- jabba(catch = catch_data, pr = priors)
#' results <- jabbaExtract(fit)
#' 
#' # List of fits (multiple scenarios)
#' fits <- list(scenario1 = fit1, scenario2 = fit2)
#' results <- jabbaExtract(fits)
#' 
#' # Nested list (multiple stocks, multiple scenarios)
#' fits <- list(
#'   stock1 = list(scenario1 = fit1a, scenario2 = fit1b),
#'   stock2 = list(scenario1 = fit2a, scenario2 = fit2b)
#' )
#' results <- jabbaExtract(fits)
#' }
#'
#' @seealso \code{\link{jabba}}, \code{\link{jabbaDiagnostics}}
#'
#' @export
setGeneric("jabbaExtract", function(object, ...) standardGeneric("jabbaExtract"))

#' @rdname jabbaExtract
#' @export
setMethod("jabbaExtract", signature(object = "ANY"),
          definition = function(object, ...) {
            if (is.null(object)) return(NULL)
            
            if (is.list(object)) {
              # Check for direct object structure first (fast path)
              if ("pars_posterior" %in% names(object)) {
                return(jabbaExtractFn(object))
              }
              if ("input" %in% names(object) && "fit" %in% names(object)) {
                return(jabbaExtractFn(object))
              }
              
              # Handle sublists (e.g., fitEB[1], fitEB[1:2])
              if (length(object) > 0) {
                # Check first element to determine structure
                firstElem <- object[[1]]
                if (!is.null(firstElem)) {
                  # Check if first element has the expected structure
                  if (is.list(firstElem)) {
                    if ("pars_posterior" %in% names(firstElem)) {
                      return(jabbaExtractList(object))
                    }
                    if (all(c("input", "fit") %in% names(firstElem))) {
                      return(jabbaExtractList(object))
                    }
                    # Check nested structure (list of lists)
                    if (length(firstElem) > 0 && !is.null(firstElem[[1]]) && 
                        is.list(firstElem[[1]]) &&
                        all(c("input", "fit") %in% names(firstElem[[1]]))) {
                      return(jabbaExtractLists(object))
                    }
                  }
                }
              }
            }
            
            return(NULL)
          })

#' Plot Summary for JABBA Fit Objects
#'
#' @description
#' Creates summary plots for a list of JABBA fit objects. This function
#' extracts the fit component from each JABBA object (which typically has
#' structure list(input, fit)) and generates diagnostic/summary plots.
#'
#' @param jb List of JABBA objects. Each element should be a list with at
#'   least two components, where the second component (\code{[[2]]}) is the
#'   JABBA fit object. Typically created by \code{jabba()} which returns
#'   \code{list(input = ..., fit = ...)}.
#'
#' @return The result of \code{jbplot_summary()} applied to the successfully
#'   extracted fit objects. Returns NULL if no valid fits are found.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Extracts the fit component (second element) from each JABBA object
#'   \item Uses \code{tryIt()} to safely handle extraction errors
#'   \item Filters out NULL results (failed extractions)
#'   \item Passes valid fits to \code{jbplot_summary()} for plotting
#' }
#'
#' The function is designed to handle lists where some elements may be NULL
#' or have invalid structure, gracefully skipping them.
#'
#' @examples
#' \dontrun{
#' # Create list of JABBA fits
#' fits <- list(
#'   scenario1 = jabba(catch = catch1, pr = priors1),
#'   scenario2 = jabba(catch = catch2, pr = priors2)
#' )
#' 
#' # Generate summary plots
#' jbplot(fits)
#' }
#'
#' @seealso \code{\link{jabba}}, \code{\link{jabbaDiagnostics}}
#'
#' @export
jbplot <- function(jb) {
  dt <- llply(jb, function(x) tryIt(x[[2]]))
  jbplot_summary(dt[!laply(dt, is.null)])
}
