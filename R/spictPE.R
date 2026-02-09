#' Extract Production Function and Process Error from SPiCT Fit
#'
#' @description
#' Extracts production function parameters and process error information from
#' a SPiCT (Stochastic Production model in Continuous Time) model fit object.
#' This function calculates both the theoretical Pella-Tomlinson production
#' curve and the realized production trajectory from the fitted model.
#'
#' @param fit A SPiCT model fit object (typically from \code{fit.spict()})
#' @param xMax Numeric. Maximum relative biomass (B/K) for production function
#'   evaluation (default: 1.3). The production function is evaluated from
#'   0.01 to \code{xMax}.
#' @param nX Numeric. Number of points at which to evaluate the production
#'   function (default: 200). Higher values provide smoother curves but
#'   increase computation time.
#'
#' @return A list with two components:
#'   \item{pf}{A data.frame containing the theoretical production function
#'     with columns:
#'     \itemize{
#'       \item \code{bRel}: Relative biomass (B/K)
#'       \item \code{b}: Absolute biomass (B)
#'       \item \code{sp}: Surplus production (P)
#'     }}
#'   \item{traj}{A data.frame containing the realized production trajectory
#'     with columns:
#'     \itemize{
#'       \item \code{year}: Year (end of interval)
#'       \item \code{sp}: Realized surplus production
#'       \item \code{spHat}: Predicted surplus production
#'       \item \code{bRel}: Starting relative biomass (B/K)
#'       \item \code{b}: Starting absolute biomass
#'       \item \code{catch}: Observed catch
#'       \item \code{cHat}: Predicted catch
#'       \item \code{bRsdl}: Biomass residuals (from fit$process.resid)
#'       \item \code{fRsdl}: Fishing mortality residuals (from fit$process.resid)
#'     }}
#'
#' @details
#' This function performs two main calculations:
#'
#' \strong{1. Production Function (Pella-Tomlinson):}
#' The theoretical production function is calculated as:
#' \deqn{P = \gamma \times MSY \times (x - x^n)}{P = gamma * MSY * (x - x^n)}
#'
#' where:
#' \itemize{
#'   \item \eqn{\gamma = \frac{1}{n-1} \times n^{n/(n-1)}}{gamma = (1/(n-1)) * n^(n/(n-1))}
#'   \item \eqn{x = B/K} (relative biomass)
#'   \item \eqn{n} is the shape parameter from SPiCT
#'   \item \eqn{MSY} is maximum sustainable yield
#'   \item \eqn{K} is carrying capacity
#' }
#'
#' \strong{2. Realized Production Trajectory:}
#' The realized production is calculated from the biomass dynamics:
#' \deqn{P_t = B_{t+1} - B_t + C_t}{P_t = B_{t+1} - B_t + C_t}
#'
#' where \eqn{C_t} is catch. This represents the actual surplus production
#' observed in the data, which can be compared to the theoretical curve.
#'
#' The function extracts:
#' \itemize{
#'   \item Model parameters (n, MSY, K) from the SPiCT fit
#'   \item Biomass estimates at integer years
#'   \item Observed and predicted catches
#'   \item Process residuals (biomass and fishing mortality)
#' }
#'
#' @note
#' This function requires the \code{spict} package and expects a SPiCT fit
#' object with the standard structure. The function uses \code{get.par()} from
#' the spict package to extract parameters.
#'
#' @examples
#' \dontrun{
#' library(spict)
#'
#' # Fit SPiCT model (example)
#' # data <- list(obsC = catch_data, timeC = catch_years, ...)
#' # fit <- fit.spict(data)
#'
#' # Extract production function and process error
#' pe_results <- spictPE(fit, xMax = 1.5, nX = 300)
#'
#' # Plot production function
#' plot(pe_results$pf$bRel, pe_results$pf$sp, type = "l",
#'      xlab = "Relative Biomass (B/K)", ylab = "Surplus Production")
#'
#' # Plot realized trajectory
#' plot(pe_results$traj$bRel, pe_results$traj$sp,
#'      xlab = "Relative Biomass", ylab = "Realized Production")
#' points(pe_results$traj$bRel, pe_results$traj$spHat, col = "red")
#' }
#'
#' @seealso
#' \code{\link[spict]{fit.spict}} for fitting SPiCT models,
#' \code{\link[spict]{get.par}} for extracting SPiCT parameters
#'
#' @export
spictPE <- function(fit, xMax = 1.3, nX = 200) {
  # Validate inputs
  if (is.null(fit)) {
    stop("'fit' cannot be NULL")
  }
  if (!is.numeric(xMax) || length(xMax) != 1 || xMax <= 0) {
    stop("'xMax' must be a positive numeric value")
  }
  if (!is.numeric(nX) || length(nX) != 1 || nX <= 0 || nX != round(nX)) {
    stop("'nX' must be a positive integer")
  }
  
  # Check if spict package is available
  if (!requireNamespace("spict", quietly = TRUE)) {
    stop("Package 'spict' is required but not installed")
  }
  
  # Extract parameters for Pella-Tomlinson curve
  # Note: get.par is from spict package
  n   <- spict::get.par("logn", fit, exp = TRUE)[1, "est"]
  msy <- spict::get.par("MSY", fit, exp = FALSE)[1, "est"]
  K   <- spict::get.par("logK", fit, exp = TRUE)[1, "est"]
  
  # Validate extracted parameters
  if (is.na(n) || is.na(msy) || is.na(K)) {
    stop("Could not extract required parameters from SPiCT fit object")
  }
  if (n <= 1) {
    stop("Shape parameter 'n' must be greater than 1 for Pella-Tomlinson model")
  }
  
  # Create relative biomass sequence
  x <- seq(0.01, xMax, length.out = nX)
  
  # Calculate gamma parameter for Pella-Tomlinson
  gamma <- (1 / (n - 1)) * n^(n / (n - 1))
  
  # Calculate production function
  Pfun <- gamma * msy * (x - x^n)
  
  # Create production function data.frame
  prod_fun <- data.frame(
    bRel = x,
    b    = x * K,
    sp   = Pfun
  )
  
  # Extract realized production trajectory
  # Biomass on SPiCT grid
  Bfull <- spict::get.par("logB", fit, exp = TRUE)[, "est"]
  tFull <- fit$inp$time
  
  # Keep only integer years
  sel <- (tFull - floor(tFull)) == 0
  B   <- Bfull[sel]
  t   <- tFull[sel]
  
  # Validate biomass extraction
  if (length(B) < 2) {
    stop("Insufficient biomass estimates for production calculation")
  }
  
  # Extract catches aligned to biomass intervals
  C <- fit$inp$obsC
  if (length(C) < length(B) - 1) {
    warning("Catch vector shorter than expected. Truncating to match biomass intervals.")
    C <- C[seq_len(min(length(C), length(B) - 1))]
  } else {
    C <- C[seq_len(length(B) - 1)]
  }
  
  # Extract predicted catches
  logCpred <- fit[[1]][names(fit[[1]]) == "logCpred"]
  if (length(logCpred) == 0) {
    warning("Predicted catches not found in fit object. Setting to NA.")
    CHat <- rep(NA, length(B) - 1)
  } else {
    CHat <- exp(logCpred)
    if (length(CHat) < length(B) - 1) {
      CHat <- CHat[seq_len(length(CHat))]
    } else {
      CHat <- CHat[seq_len(length(B) - 1)]
    }
  }
  
  # Calculate realized production: P_t = B_{t+1} - B_t + C_t
  P    <- B[-1] - B[-length(B)] + C
  PHat <- B[-1] - B[-length(B)] + CHat
  
  # Create production trajectory data.frame
  prod_traj <- data.frame(
    year = utils::tail(t, -1),          # end of interval
    sp   = P,
    spHat = PHat,
    bRel = B[-length(B)] / K,            # starting biomass of interval, relative to K
    b    = B[-length(B)],
    catch = C,
    cHat = CHat
  )
  
  # Add process residuals if available
  if (!is.null(fit$process.resid) && is.data.frame(fit$process.resid)) {
    # Ensure process.resid has correct column names
    if (ncol(fit$process.resid) >= 3) {
      names(fit$process.resid)[1:3] <- c("year", "bRsdl", "fRsdl")
      # Merge with production trajectory
      prod_traj <- merge(prod_traj, fit$process.resid, by = "year", all.x = TRUE)
    } else {
      warning("Process residuals structure unexpected. Skipping merge.")
    }
  } else {
    warning("Process residuals not found in fit object. Trajectory will not include residuals.")
  }
  
  # Return list with both components
  return(list(pf = prod_fun, traj = prod_traj))
}
