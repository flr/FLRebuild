#' Extract Production Function and Process Error from JABBA Fit
#'
#' @description
#' Extracts production function parameters and process error information from
#' a JABBA (Just Another Bayesian Biomass Assessment) model fit object.
#' This function calculates both the theoretical Pella-Tomlinson production
#' curve and the realized production trajectory from the fitted model.
#'
#' @param fit A JABBA model fit object (typically from \code{jabba()} or
#'   \code{fit_jabba()})
#' @param xMax Numeric. Maximum relative biomass (B/K) for production function
#'   evaluation (default: 1.5). The production function is evaluated from
#'   0.01 to \code{xMax}.
#' @param nX Numeric. Number of points at which to evaluate the production
#'   function (default: 200). Higher values provide smoother curves but
#'   increase computation time.
#' @param use Character. Which posterior statistic to use for parameters:
#'   "median" (default) or "mean"
#' @param m.default Numeric. Default value for shape parameter m if not found
#'   in posterior (default: 2, corresponding to Schaefer model)
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
#'       \item \code{bRel}: Starting relative biomass (B/K)
#'       \item \code{b}: Starting absolute biomass
#'       \item \code{catch}: Observed catch
#'       \item \code{pe}: Process error (log residuals from state equation)
#'     }}
#'
#' @details
#' This function performs two main calculations:
#'
#' \strong{1. Production Function (Pella-Tomlinson):}
#' The theoretical production function is calculated as:
#' \deqn{P = r \times B \times (1 - (B/K)^m)}{P = r * B * (1 - (B/K)^m)}
#'
#' where:
#' \itemize{
#'   \item \eqn{r} is the intrinsic growth rate from JABBA
#'   \item \eqn{B} is absolute biomass
#'   \item \eqn{K} is carrying capacity from JABBA
#'   \item \eqn{m} is the shape parameter (m=2 for Schaefer, m=1 for Fox)
#'   \item \eqn{x = B/K} (relative biomass)
#' }
#'
#' \strong{2. Realized Production Trajectory:}
#' The realized production is calculated from the biomass dynamics:
#' \deqn{P_t = B_{t+1} - B_t + C_t}{P_t = B_{t+1} - B_t + C_t}
#'
#' where \eqn{C_t} is catch. This represents the actual surplus production
#' observed in the data, which can be compared to the theoretical curve.
#'
#' Process error is calculated as the log residual from the state equation:
#' \deqn{\epsilon_t = \log(B_{t+1} / \mu_t)}{eps_t = log(B_{t+1} / mu_t)}
#'
#' where \eqn{\mu_t = B_t + P_t - C_t}{mu_t = B_t + P_t - C_t} is the expected
#' biomass at time t+1.
#'
#' The function extracts:
#' \itemize{
#'   \item Model parameters (r, K, m) from JABBA posterior distributions
#'   \item Biomass estimates from timeseries
#'   \item Observed catches
#'   \item Process error residuals
#' }
#'
#' @note
#' This function requires a JABBA fit object with the standard structure:
#' \itemize{
#'   \item \code{fit$pars_posterior}: Matrix of posterior parameter samples
#'   \item \code{fit$timeseries}: Array with biomass timeseries
#'   \item \code{fit$est.catch}: Matrix with estimated catches
#' }
#'
#' @examples
#' \dontrun{
#' library(JABBA)
#'
#' # Fit JABBA model (example)
#' # jb <- jabba(data, ...)
#'
#' # Extract production function and process error
#' pe_results <- jabbaPE(jb, xMax = 1.5, nX = 300)
#'
#' # Plot production function
#' plot(pe_results$pf$bRel, pe_results$pf$sp, type = "l",
#'      xlab = "Relative Biomass (B/K)", ylab = "Surplus Production")
#'
#' # Plot realized trajectory
#' plot(pe_results$traj$bRel, pe_results$traj$sp,
#'      xlab = "Relative Biomass", ylab = "Realized Production")
#' }
#'
#' @seealso
#' \code{\link{spictPE}} for similar function for SPiCT models
#'
#' @export
jabbaPE <- function(fit, xMax = 1.5, nX = 200, use = c("median", "mean"), m.default = 2) {
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
  use <- match.arg(use)
  
  # Check required JABBA structure
  if (is.null(fit$pars_posterior)) {
    stop("JABBA fit object must have 'pars_posterior' component")
  }
  if (is.null(fit$timeseries)) {
    stop("JABBA fit object must have 'timeseries' component")
  }
  if (is.null(fit$est.catch)) {
    stop("JABBA fit object must have 'est.catch' component")
  }
  
  # Extract parameters from posterior
  parpost <- fit$pars_posterior
  
  get_par <- function(name) {
    if (!name %in% colnames(parpost)) {
      if (name == "m") {
        return(m.default)
      }
      stop("Parameter '", name, "' not found in fit$pars_posterior")
    }
    if (use == "median") {
      median(parpost[, name])
    } else {
      mean(parpost[, name])
    }
  }
  
  # Extract parameters
  r <- get_par("r")
  K <- get_par("K")
  m <- get_par("m")
  
  # Validate parameters
  if (is.na(r) || is.na(K) || is.na(m)) {
    stop("Could not extract required parameters (r, K, m) from JABBA fit object")
  }
  if (r <= 0 || K <= 0 || m <= 0) {
    stop("Parameters r, K, and m must be positive")
  }
  
  # Extract timeseries data
  ts <- fit$timeseries
  B <- ts[, 1, "B"]
  C <- fit$est.catch[, "mu"]
  years <- as.numeric(dimnames(ts)[[1]])
  
  # Validate data
  n <- length(B)
  if (length(C) != n) {
    stop("Biomass and catch vectors have different lengths. Check fit$timeseries and fit$est.catch structure.")
  }
  if (n < 2) {
    stop("Insufficient time series data for production calculation")
  }
  
  # Create relative biomass sequence for production function
  x <- seq(0.01, xMax, length.out = nX)
  
  # Calculate theoretical production function: P = r * B * (1 - (B/K)^m)
  # where B = x * K, so: P = r * x * K * (1 - x^m)
  Pfun <- r * x * K * (1 - x^m)
  
  # Create production function data.frame
  prod_fun <- data.frame(
    bRel = x,
    b    = x * K,
    sp   = Pfun
  )
  
  # Calculate realized production trajectory: P_t = B_{t+1} - B_t + C_t
  P <- B[-1] - B[-length(B)] + C[-length(C)]
  
  # Calculate process error: eps_t = log(B_{t+1} / mu_t)
  # where mu_t = B_t + P_t - C_t
  mu <- rep(NA_real_, n - 1)
  eps <- rep(NA_real_, n - 1)
  
  for (t in seq_len(n - 1)) {
    Bt <- B[t]
    Ct <- C[t]
    
    # Expected production
    prod <- r * Bt * (1 - (Bt / K)^m)
    
    # Expected next biomass
    mu[t] <- Bt + prod - Ct
    mu[t] <- max(mu[t], 1e-12)  # Ensure positive
    
    # Process error (log residual)
    eps[t] <- log(B[t + 1] / mu[t])
  }
  
  # Create production trajectory data.frame
  prod_traj <- data.frame(
    year  = years[-1],              # end of interval
    sp    = P,
    bRel  = B[-length(B)] / K,      # starting biomass of interval, relative to K
    b     = B[-length(B)],
    catch = C[-length(C)],
    pe    = eps
  )
  
  # Return list with both components (consistent with spictPE structure)
  return(list(pf = prod_fun, traj = prod_traj))
}
