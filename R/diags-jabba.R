# Constants
P_VALUE_THRESHOLD  = 0.05
BIAS_THRESHOLD     = 0.1
SHAPIRO_MIN_SAMPLE = 3
SHAPIRO_MAX_SAMPLE = 5000
LJUNG_BOX_LAG      = 5

#' Extract P-Values from Matrix or Vector
#'
#' @description
#' Helper function to extract p-values from convergence test results.
#' If input is a matrix, extracts the first column; otherwise returns
#' the vector as-is.
#'
#' @param pVals Numeric matrix or vector. P-values from convergence tests
#'   (e.g., Geweke, Heidelberger-Welch)
#'
#' @return Numeric vector of p-values
#'
#' @keywords internal
getPValues<-function(pVals) {
  if (is.matrix(pVals)) {
    pVals[, 1]
  } else {
    pVals}}

#' Check Convergence Test Results
#'
#' @description
#' Helper function to evaluate convergence test results. Determines if all
#' p-values exceed the threshold (test passed) and calculates the mean
#' p-value across all parameters.
#'
#' @param pVals Numeric matrix or vector. P-values from convergence tests
#'   (e.g., Geweke, Heidelberger-Welch)
#' @param threshold Numeric. P-value threshold for passing test. Default:
#'   \code{P_VALUE_THRESHOLD} (0.05)
#'
#' @return List with components:
#'   \itemize{
#'     \item \code{passed}: Logical. TRUE if all p-values > threshold
#'     \item \code{mean}: Numeric. Mean p-value across all parameters
#'   }
#'   Returns \code{list(passed = NA, mean = NA)} if \code{pVals} is NULL.
#'
#' @keywords internal
checkConvergenceTest<-function(pVals, threshold = P_VALUE_THRESHOLD) {
  if (is.null(pVals)) {
    return(list(passed = NA, mean = NA))
  }
  vals = getPValues(pVals)
  list(
    passed = all(vals > threshold, na.rm = TRUE),
    mean = mean(vals, na.rm = TRUE))}

#' Get Convergence Diagnostics from JABBA Fit
#'
#' @description
#' Extracts and evaluates convergence diagnostics from a JABBA fit object.
#' Performs Geweke and Heidelberger-Welch convergence tests on MCMC chains.
#'
#' @param jb JABBA fit object. Must have \code{pars} component containing
#'   \code{Geweke.p} and \code{Heidel.p} p-values.
#'
#' @return List with components:
#'   \itemize{
#'     \item \code{Geweke_passed}: Logical. Whether Geweke test passed
#'       (all p-values > 0.05)
#'     \item \code{Geweke_mean}: Numeric. Mean Geweke p-value
#'     \item \code{Heidel_passed}: Logical. Whether Heidelberger-Welch test passed
#'       (all p-values > 0.05)
#'     \item \code{Heidel_mean}: Numeric. Mean Heidel p-value
#'   }
#'   Returns NAs if \code{jb$pars} is NULL.
#'
#' @details
#' Convergence diagnostics test whether MCMC chains have converged to the
#' target distribution:
#' \itemize{
#'   \item \strong{Geweke test:} Tests for convergence by comparing means
#'     from the beginning and end of chains
#'   \item \strong{Heidelberger-Welch test:} Tests for stationarity of chains
#' }
#' Both tests should have p-values > 0.05 for convergence.
#'
#' @seealso \code{\link{getResidualDiagnostics}}, \code{\link{getTerminalParameters}}
#'
#' @keywords internal
getConvergenceDiagnostics<-function(jb,pars=c("K","r","m")) {
  result = list(
    Geweke_passed = NA,
    Geweke_mean   = NA,
    Heidel_passed = NA,
    Heidel_mean   = NA)
  
  if (!is.null(jb$pars)) {
    geweke = checkConvergenceTest(jb$pars[pars,"Geweke.p"])
    result$Geweke_passed = geweke$passed
    result$Geweke_mean   = geweke$mean
    
    heidel = checkConvergenceTest(jb$pars[pars,"Heidel.p"])
    result$Heidel_passed = heidel$passed
    result$Heidel_mean   = heidel$mean
  }
  
  return(result)}

#' Get Residual Diagnostics from JABBA Fit
#'
#' @description
#' Extracts and evaluates residual diagnostics from a JABBA fit object.
#' Performs statistical tests on model residuals including normality,
#' autocorrelation, runs test, and bias checks.
#'
#' @param jb JABBA fit object. Must have \code{diag} component containing
#'   \code{residual} vector.
#'
#' @return List with components:
#'   \itemize{
#'     \item \code{n_residuals}: Integer. Number of residuals
#'     \item \code{mean_residual}: Numeric. Mean residual value
#'     \item \code{sd_residual}: Numeric. Standard deviation of residuals
#'     \item \code{Normality_passed}: Logical. Whether residuals pass normality
#'       test (Shapiro-Wilk, p > 0.05). Only calculated if 3 <= n <= 5000
#'     \item \code{Normality_p}: Numeric. Shapiro-Wilk p-value
#'     \item \code{Autocorr_passed}: Logical. Whether residuals pass
#'       autocorrelation test (Ljung-Box, p > 0.05)
#'     \item \code{Autocorr_p}: Numeric. Ljung-Box p-value (lag = 5)
#'     \item \code{Runs_passed}: Logical. Whether runs test passed (p > 0.05)
#'     \item \code{Runs_p}: Numeric. Runs test p-value
#'     \item \code{Bias_passed}: Logical. Whether mean residual bias is
#'       acceptable (|mean| < 0.1)
#'   }
#'   Returns NAs if \code{jb$diag$residual} is NULL or missing.
#'
#' @details
#' Residual diagnostics test model assumptions and fit quality:
#' \itemize{
#'   \item \strong{Normality:} Shapiro-Wilk test (requires 3-5000 residuals)
#'   \item \strong{Autocorrelation:} Ljung-Box test with lag 5
#'   \item \strong{Runs test:} Tests for randomness in residual sequence
#'     (requires \code{jbrunstest} function)
#'   \item \strong{Bias:} Checks if mean residual is close to zero
#' }
#'
#' @seealso \code{\link{getConvergenceDiagnostics}}, \code{\link{getTerminalParameters}},
#'   \code{\link[stats]{shapiro.test}}, \code{\link[stats]{Box.test}}
#'
#' @keywords internal
getResidualDiagnostics<-function(jb) {
  result = list(
    n_residuals      = NA,
    mean_residual    = NA,
    sd_residual      = NA,
    Normality_passed = NA,
    Normality_p      = NA,
    Autocorr_passed  = NA,
    Autocorr_p       = NA,
    Runs_passed      = NA,
    Runs_p           = NA,
    Bias_passed      = NA)
  
  if (is.null(jb$diag) || !"residual" %in% names(jb$diag)) {
    return(result)}
  
  residuals  = jb$diag$residual
  nResiduals = length(residuals)
  
  result$n_residuals   = nResiduals
  result$mean_residual = mean(residuals, na.rm = TRUE)
  result$sd_residual   = sd(residuals,   na.rm = TRUE)
  
  # Normality test (Shapiro-Wilk)
  if (nResiduals >= SHAPIRO_MIN_SAMPLE && nResiduals <= SHAPIRO_MAX_SAMPLE) {
    shapiroTest = shapiro.test(residuals)
    result$Normality_passed = shapiroTest$p.value > P_VALUE_THRESHOLD
    result$Normality_p = shapiroTest$p.value}
  
  # Autocorrelation test (Ljung-Box)
  ljungBox = Box.test(residuals, lag = LJUNG_BOX_LAG, type = "Ljung-Box")
  result$Autocorr_passed = ljungBox$p.value > P_VALUE_THRESHOLD
  result$Autocorr_p = ljungBox$p.value
  
  # Runs test
  runsTest = jbrunstest(jb)
  result$Runs_passed = runsTest$runs.p > P_VALUE_THRESHOLD
  result$Runs_p = runsTest$runs.p
  
  # Bias test
  result$Bias_passed = abs(result$mean_residual) < BIAS_THRESHOLD
  
  return(result)}

#' Get Terminal Parameter Estimates from JABBA Fit
#'
#' @description
#' Extracts terminal (final year) stock status indicators from a JABBA fit
#' object. Retrieves B/BMSY and F/FMSY ratios from the Kobe matrix.
#'
#' @param jb JABBA fit object. Must have \code{kobe} component containing
#'   a data.frame with \code{stock} (B/BMSY) and/or \code{harvest} (F/FMSY)
#'   columns.
#'
#' @return List with components:
#'   \itemize{
#'     \item \code{BBmsy_terminal}: Numeric. Terminal B/BMSY ratio (relative
#'       stock biomass)
#'     \item \code{FFmsy_terminal}: Numeric. Terminal F/FMSY ratio (relative
#'       fishing mortality)
#'   }
#'   Returns NAs if \code{jb$kobe} is NULL, not a data.frame, or empty.
#'
#' @details
#' Terminal parameters indicate the final year stock status:
#' \itemize{
#'   \item \strong{B/BMSY:} Values < 1 indicate overfished stock
#'   \item \strong{F/FMSY:} Values > 1 indicate overfishing
#' }
#' These values are extracted from the last row of the Kobe matrix, which
#' contains the posterior distributions of stock and harvest status over time.
#'
#' @seealso \code{\link{getConvergenceDiagnostics}}, \code{\link{getResidualDiagnostics}}
#'
#' @keywords internal
getTerminalParameters<-function(jb) {
  result = list(
    BBmsy_terminal = NA,
    FFmsy_terminal = NA)
  
  if (is.null(jb$kobe) || !is.data.frame(jb$kobe) || nrow(jb$kobe) == 0) {
    return(result)}
  
  terminalKobe = jb$kobe[nrow(jb$kobe), ]
  
  if ("stock" %in% names(terminalKobe)) {
    result$BBmsy_terminal = terminalKobe$stock[1]}
  if ("harvest" %in% names(terminalKobe)) {
    result$FFmsy_terminal = terminalKobe$harvest[1]}
  
  return(result)}

#' Create Empty Diagnostic Result Data Frame
#'
#' @description
#' Helper function to create an empty diagnostic result data.frame with
#' all required columns initialized to NA. Used when diagnostics cannot be
#' calculated (e.g., missing data, fit errors).
#'
#' @param id Character string. Stock identifier
#' @param scenario Character string. Scenario identifier
#'
#' @return data.frame with columns:
#'   \itemize{
#'     \item \code{Stock}: Stock identifier
#'     \item \code{Scenario}: Scenario identifier
#'     \item All diagnostic columns initialized to NA (Geweke_passed, Geweke_mean,
#'       Heidel_passed, Heidel_mean, n_residuals, mean_residual, sd_residual,
#'       Normality_passed, Normality_p, Autocorr_passed, Autocorr_p,
#'       Runs_passed, Runs_p, Bias_passed, BBmsy_terminal, FFmsy_terminal)
#'   }
#'
#' @keywords internal
createEmptyResult<-function(id, scenario) {
  data.frame(
    Stock = id,
    Scenario = scenario,
    Geweke_passed    = NA,
    Geweke_mean      = NA,
    Heidel_passed    = NA,
    Heidel_mean      = NA,
    n_residuals      = NA,
    mean_residual    = NA,
    sd_residual      = NA,
    Normality_passed = NA,
    Normality_p      = NA,
    Autocorr_passed  = NA,
    Autocorr_p       = NA,
    Runs_passed      = NA,
    Runs_p           = NA,
    Bias_passed      = NA,
    BBmsy_terminal   = NA,
    FFmsy_terminal   = NA,
    stringsAsFactors = FALSE)}

#' Run All Diagnostics for a Single JABBA Fit
#'
#' @description
#' Runs comprehensive diagnostics on a single JABBA fit object. Combines
#' convergence diagnostics, residual diagnostics, and terminal parameter
#' estimates into a single data.frame.
#'
#' @param jb JABBA fit object. Should have components: \code{pars} (for
#'   convergence), \code{diag} (for residuals), and \code{kobe} (for terminal
#'   parameters).
#'
#' @return data.frame with all diagnostic results combined:
#'   \itemize{
#'     \item Convergence diagnostics (Geweke, Heidel)
#'     \item Residual diagnostics (normality, autocorrelation, runs, bias)
#'     \item Terminal parameters (B/BMSY, F/FMSY)
#'   }
#'   Missing components result in NA values for corresponding diagnostics.
#'
#' @details
#' This function is the main entry point for running diagnostics on a single
#' JABBA fit. It calls three diagnostic functions and combines their results:
#' \itemize{
#'   \item \code{getConvergenceDiagnostics()}: MCMC convergence tests
#'   \item \code{getResidualDiagnostics()}: Residual statistical tests
#'   \item \code{getTerminalParameters()}: Final year stock status
#' }
#'
#' @seealso \code{\link{runAll}}, \code{\link{getConvergenceDiagnostics}},
#'   \code{\link{getResidualDiagnostics}}, \code{\link{getTerminalParameters}}
#'
#' @keywords internal
runDiags<-function(jb,pars=c("K","r","m")) {
    
  # get all diagnostics
  convergence = getConvergenceDiagnostics(jb,pars=pars)
  residuals   = getResidualDiagnostics(jb)
  terminal    = getTerminalParameters(jb)
    
  # Combine all results
  result = cbind(as.data.frame(convergence), 
                 as.data.frame(residuals), 
                 as.data.frame(terminal))
    
  return(result)}


#' Standardize Columns Across Diagnostic Data Frames
#'
#' @description
#' Helper function to ensure all diagnostic data.frames in a list have the
#' same columns. Missing columns are added and filled with NA values. This is
#' necessary when combining diagnostics from multiple fits that may have
#' different available diagnostics.
#'
#' @param diagsList List of data.frames. Each element should be a diagnostic
#'   result data.frame (typically from \code{runDiags()}).
#'
#' @return List of data.frames with identical column structure. All data.frames
#'   will have all columns that appear in any element, with missing columns
#'   filled with NA.
#'
#' @details
#' When combining diagnostics from multiple JABBA fits, some fits may have
#' missing components (e.g., no residuals, no convergence diagnostics).
#' This function ensures all data.frames can be safely combined using
#' \code{rbind()} by standardizing their column structure.
#'
#' @keywords internal
chkCols<-function(diagsList) {
  allCols = unique(unlist(lapply(diagsList, names)))
   
  lapply(diagsList, function(x) {
    missingCols = setdiff(allCols, names(x))
    if (length(missingCols) > 0) {
      for (col in missingCols) {
        x[[col]] = NA}}
    x[, allCols, drop = FALSE]})}

#' Calculate Summary Statistics for Diagnostic Results
#'
#' @description
#' Calculates summary statistics across multiple diagnostic results (e.g.,
#' for a scenario with multiple stocks). Counts how many stocks passed each
#' diagnostic test and total number of tests performed.
#'
#' @param data data.frame. Diagnostic results with columns for each diagnostic
#'   test (Geweke_passed, Heidel_passed, Normality_passed, etc.). Typically
#'   from \code{runAll()} or combined results from multiple \code{runDiags()}
#'   calls.
#'
#' @return data.frame with summary statistics:
#'   \itemize{
#'     \item \code{nStocks}: Total number of stocks
#'     \item \code{nGewekePassed}: Number passing Geweke test
#'     \item \code{nHeidelPassed}: Number passing Heidel test
#'     \item \code{nNormalPassed}: Number passing normality test
#'     \item \code{nNormalTotal}: Total number with normality test results
#'     \item \code{nAutocorrPassed}: Number passing autocorrelation test
#'     \item \code{nAutocorrTotal}: Total number with autocorrelation test results
#'     \item \code{nRunsPassed}: Number passing runs test
#'     \item \code{nRunsTotal}: Total number with runs test results
#'     \item \code{nBiasPassed}: Number passing bias test
#'     \item \code{nBiasTotal}: Total number with bias test results
#'   }
#'
#' @details
#' This function aggregates diagnostic results across multiple stocks to
#' provide scenario-level summary statistics. It counts both the number of
#' passing tests and the total number of tests performed (since some tests
#' may not be calculable for all stocks).
#'
#' @seealso \code{\link{runAll}}, \code{\link{runDiags}}
#'
#' @keywords internal
calcSmry<-function(data) {
  data.frame(
    nStocks = nrow(data),
    nGewekePassed   = sum(data$Geweke_passed, na.rm = TRUE),
    nHeidelPassed   = sum(data$Heidel_passed, na.rm = TRUE),
    nNormalPassed   = sum(data$Normality_passed, na.rm = TRUE),
    nNormalTotal    = sum(!is.na(data$Normality_passed)),
    nAutocorrPassed = sum(data$Autocorr_passed, na.rm = TRUE),
    nAutocorrTotal  = sum(!is.na(data$Autocorr_passed)),
    nRunsPassed     = sum(data$Runs_passed, na.rm = TRUE),
    nRunsTotal      = sum(!is.na(data$Runs_passed)),
    nBiasPassed     = sum(data$Bias_passed, na.rm = TRUE),
    nBiasTotal      = sum(!is.na(data$Bias_passed)))}

#' Run Diagnostics for Multiple JABBA Fits
#'
#' @description
#' Runs diagnostics on a list of JABBA fit objects (e.g., multiple stocks
#' or scenarios). Processes each fit and combines results into a single
#' data.frame with standardized columns.
#'
#' @param jb Named list of JABBA fit objects. Each element should be a JABBA
#'   fit object that can be processed by \code{runDiags()}. List names are
#'   used as identifiers (typically stock IDs or scenario names).
#'
#' @return data.frame with diagnostic results for all fits, combined using
#'   \code{rbind()}. Contains all diagnostic columns (convergence, residuals,
#'   terminal parameters). Fits that return NULL are skipped. All data.frames
#'   are standardized to have the same columns before combining.
#'
#' @details
#' This function processes multiple JABBA fits in a single call:
#' \itemize{
#'   \item Loops through each element in the list
#'   \item Calls \code{runDiags()} for each fit
#'   \item Skips NULL results (failed diagnostics)
#'   \item Standardizes columns using \code{chkCols()} to ensure all
#'     data.frames can be combined
#'   \item Combines all results into a single data.frame
#' }
#'
#' This is useful for batch processing diagnostics across multiple stocks
#' or scenarios.
#'
#' @examples
#' \dontrun{
#' # Run diagnostics for multiple stocks
#' jb_fits <- list(
#'   stock1 = jabba_fit1,
#'   stock2 = jabba_fit2,
#'   stock3 = jabba_fit3
#' )
#' all_diagnostics <- runAll(jb_fits)
#' }
#'
#' @seealso \code{\link{runDiags}}, \code{\link{chkCols}}, \code{\link{calcSmry}}
#'
#' @keywords internal
runAll<-function(jb,pars=c("K","r","m")){
  rtn = list()
  for (id in names(jb)) {
    rslt = runDiags(jb[[id]],pars=pars)
    if (!is.null(rslt))
      rtn[[id]] = rslt}
  
  # Standardize columns after each stock
  do.call(rbind, chkCols(rtn))}

# jbICES=llply(histICES, function(x) x[["ICES"]]$fit)
# dgICES=runAll(jbICES)
