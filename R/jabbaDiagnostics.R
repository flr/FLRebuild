#' Run Diagnostics on JABBA Model Fits
#'
#' @description
#' Generic function to run comprehensive diagnostics on JABBA (Just Another
#' Bayesian Biomass Assessment) model fits. This includes convergence tests,
#' residual diagnostics, and parameter estimates. The function is designed to
#' work with single fits or lists of fits.
#'
#' @param fit A JABBA fit object (list with $fit element containing the JABBA model)
#' @param id Character. Identifier for the stock/model (optional, for tracking)
#' @param scenario Character. Scenario identifier (optional, for tracking)
#' @param ... Additional arguments (currently unused)
#'
#' @return A data.frame with diagnostic results containing:
#'   \item{Stock}{Stock identifier (if provided)}
#'   \item{Scenario}{Scenario identifier (if provided)}
#'   \item{Geweke_passed}{Logical. Whether Geweke convergence test passed (all p > 0.05)}
#'   \item{Geweke_mean}{Numeric. Mean Geweke p-value}
#'   \item{Heidel_passed}{Logical. Whether Heidelberger-Welch test passed (all p > 0.05)}
#'   \item{Heidel_mean}{Numeric. Mean Heidel p-value}
#'   \item{n_residuals}{Integer. Number of residuals}
#'   \item{mean_residual}{Numeric. Mean residual value}
#'   \item{sd_residual}{Numeric. Standard deviation of residuals}
#'   \item{Normality_passed}{Logical. Whether residuals pass normality test (Shapiro-Wilk, p > 0.05)}
#'   \item{Normality_p}{Numeric. Shapiro-Wilk p-value}
#'   \item{Autocorr_passed}{Logical. Whether residuals pass autocorrelation test (Ljung-Box, p > 0.05)}
#'   \item{Autocorr_p}{Numeric. Ljung-Box p-value}
#'   \item{Runs_passed}{Logical. Whether runs test passed (p > 0.05)}
#'   \item{Runs_p}{Numeric. Runs test p-value}
#'   \item{Bias_passed}{Logical. Whether mean residual bias is acceptable (|mean| < 0.1)}
#'   \item{BBmsy_terminal}{Numeric. Terminal B/Bmsy ratio}
#'   \item{FFmsy_terminal}{Numeric. Terminal F/Fmsy ratio}
#'
#' @details
#' This function performs three categories of diagnostics:
#'
#' \strong{1. Convergence Diagnostics:}
#' \itemize{
#'   \item \strong{Geweke test:} Tests for convergence of MCMC chains
#'   \item \strong{Heidelberger-Welch test:} Tests for stationarity of MCMC chains
#' }
#'
#' \strong{2. Residual Diagnostics:}
#' \itemize{
#'   \item \strong{Normality:} Shapiro-Wilk test for residual normality (n between 3-5000)
#'   \item \strong{Autocorrelation:} Ljung-Box test for residual autocorrelation (lag=5)
#'   \item \strong{Runs test:} Tests for randomness in residual sequence (if jbrunstest available)
#'   \item \strong{Bias:} Checks if mean residual is close to zero (|mean| < 0.1)
#' }
#'
#' \strong{3. Parameter Estimates:}
#' \itemize{
#'   \item Terminal B/Bmsy and F/Fmsy from kobe matrix
#' }
#'
#' The function handles missing data gracefully, returning NA for tests that
#' cannot be performed.
#'
#' @examples
#' \dontrun{
#' # Single JABBA fit
#' diagnostics <- jabbaDiagnostics(jabba_fit, id = "cod.27.47d20", scenario = "EBiomass")
#'
#' # List of fits (nested structure: fits[[stock]][[scenario]])
#' fits <- list(
#'   stock1 = list(scenario1 = list(fit = jabba_fit1)),
#'   stock2 = list(scenario2 = list(fit = jabba_fit2))
#' )
#' all_diagnostics <- jabbaDiagnostics(fits)
#' }
#'
#' @seealso
#' \code{\link{shapiro.test}}, \code{\link{Box.test}} for statistical tests,
#' \code{\link{summarizeDiagnostics}} for summarizing results across fits
#'
#' @importFrom stats shapiro.test Box.test
#' @export
setGeneric("jabbaDiagnostics", function(fit, id = NULL, scenario = NULL, ...) {
  standardGeneric("jabbaDiagnostics")
})

#' @rdname jabbaDiagnostics
#' @export
setMethod("jabbaDiagnostics", signature(fit = "list"),
  function(fit, id = NULL, scenario = NULL, ...) {
    
    # Handle nested list structure: fits[[id]][[scenario]]$fit
    # Check structure by examining first element
    if (length(fit) == 0) {
      return(NULL)
    }
    
    first_elem <- fit[[1]]
    
    # Case 1: Multiple stocks with scenarios: fit[[stock]][[scenario]]$fit
    # First element is a list, and its first element has $fit
    if (is.list(first_elem) && length(first_elem) > 0) {
      first_scenario <- first_elem[[1]]
      if (is.list(first_scenario) && !is.null(first_scenario$fit)) {
        # This is multiple stocks: fit[[stock]][[scenario]]$fit
        results <- list()
        for (stock_id in names(fit)) {
          stock_fits <- fit[[stock_id]]
          if (is.list(stock_fits)) {
            for (scen in names(stock_fits)) {
              if (!is.null(stock_fits[[scen]]) && 
                  is.list(stock_fits[[scen]]) && 
                  !is.null(stock_fits[[scen]]$fit)) {
                result <- jabbaDiagnostics(stock_fits[[scen]]$fit, 
                                          id = stock_id, 
                                          scenario = scen, ...)
                if (!is.null(result)) {
                  results[[paste(stock_id, scen, sep = "_")]] <- result
                }
              }
            }
          }
        }
        if (length(results) > 0) {
          # Ensure all results have same columns
          all_cols <- unique(unlist(lapply(results, names)))
          results <- lapply(results, function(x) {
            missing_cols <- setdiff(all_cols, names(x))
            if (length(missing_cols) > 0) {
              for (col in missing_cols) {
                x[[col]] <- NA
              }
            }
            x[, all_cols, drop = FALSE]
          })
          return(do.call(rbind, results))
        } else {
          return(NULL)
        }
      }
    }
    
    # Case 2: Single stock with multiple scenarios: fit[[scenario]]$fit
    # First element has $fit directly
    if (is.list(first_elem) && !is.null(first_elem$fit)) {
      results <- list()
      for (scen in names(fit)) {
        if (!is.null(fit[[scen]]) && 
            is.list(fit[[scen]]) && 
            !is.null(fit[[scen]]$fit)) {
          result <- jabbaDiagnostics(fit[[scen]]$fit, id = id, scenario = scen, ...)
          if (!is.null(result)) {
            results[[scen]] <- result
          }
        }
      }
      if (length(results) > 0) {
        return(do.call(rbind, results))
      } else {
        return(NULL)
      }
    }
    
    # Case 3: Single fit object - extract $fit if present
    if (!is.null(fit$fit)) {
      jb <- fit$fit
    } else {
      jb <- fit
    }
    return(jabbaDiagnosticsSingle(jb, id = id, scenario = scenario, ...))
  })

#' Internal function to run diagnostics on a single JABBA fit object
#'
#' @param jb A JABBA model object (the $fit element)
#' @param id Character. Stock identifier
#' @param scenario Character. Scenario identifier
#' @return data.frame with diagnostic results
#' @noRd
jabbaDiagnosticsSingle <- function(jb, id = NULL, scenario = NULL) {
  
  tryCatch({
    # Validate input
    if (is.null(jb) || !is.list(jb)) {
      return(NULL)
    }
    
    # Initialize result data.frame
    result <- data.frame(
      Stock = if (is.null(id)) NA_character_ else id,
      Scenario = if (is.null(scenario)) NA_character_ else scenario,
      stringsAsFactors = FALSE
    )
    
    # 1. Convergence diagnostics
    if (!is.null(jb$pars)) {
      # Geweke test
      if (!is.null(jb$pars$Geweke.p)) {
        geweke_vals <- if (is.matrix(jb$pars$Geweke.p)) {
          jb$pars$Geweke.p[, 1]
        } else {
          jb$pars$Geweke.p
        }
        result$Geweke_passed <- all(geweke_vals > 0.05, na.rm = TRUE)
        result$Geweke_mean <- mean(geweke_vals, na.rm = TRUE)
      } else {
        result$Geweke_passed <- NA
        result$Geweke_mean <- NA
      }
      
      # Heidelberger-Welch test
      if (!is.null(jb$pars$Heidel.p)) {
        heidel_vals <- if (is.matrix(jb$pars$Heidel.p)) {
          jb$pars$Heidel.p[, 1]
        } else {
          jb$pars$Heidel.p
        }
        result$Heidel_passed <- all(heidel_vals > 0.05, na.rm = TRUE)
        result$Heidel_mean <- mean(heidel_vals, na.rm = TRUE)
      } else {
        result$Heidel_passed <- NA
        result$Heidel_mean <- NA
      }
    } else {
      result$Geweke_passed <- NA
      result$Geweke_mean <- NA
      result$Heidel_passed <- NA
      result$Heidel_mean <- NA
    }
    
    # 2. Residual diagnostics
    if (!is.null(jb$diag) && "residual" %in% names(jb$diag)) {
      residuals <- jb$diag$residual
      n_residuals <- length(residuals)
      
      result$n_residuals <- n_residuals
      result$mean_residual <- mean(residuals, na.rm = TRUE)
      result$sd_residual <- sd(residuals, na.rm = TRUE)
      
      # Normality test (Shapiro-Wilk)
      if (n_residuals >= 3 && n_residuals <= 5000) {
        shapiro_test <- shapiro.test(residuals)
        result$Normality_passed <- shapiro_test$p.value > 0.05
        result$Normality_p <- shapiro_test$p.value
      } else {
        result$Normality_passed <- NA
        result$Normality_p <- NA
      }
      
      # Autocorrelation test (Ljung-Box, lag=5)
      ljung_box <- Box.test(residuals, lag = 5, type = "Ljung-Box")
      result$Autocorr_passed <- ljung_box$p.value > 0.05
      result$Autocorr_p <- ljung_box$p.value
      
      # Runs test (if available)
      if (exists("jbrunstest", mode = "function")) {
        tryCatch({
          runs_test <- jbrunstest(jb)
          result$Runs_passed <- runs_test$runs.p > 0.05
          result$Runs_p <- runs_test$runs.p
        }, error = function(e) {
          result$Runs_passed <- NA
          result$Runs_p <- NA
        })
      } else {
        result$Runs_passed <- NA
        result$Runs_p <- NA
      }
      
      # Bias test
      result$Bias_passed <- abs(mean(residuals, na.rm = TRUE)) < 0.1
      
    } else {
      result$n_residuals <- NA
      result$mean_residual <- NA
      result$sd_residual <- NA
      result$Normality_passed <- NA
      result$Normality_p <- NA
      result$Autocorr_passed <- NA
      result$Autocorr_p <- NA
      result$Runs_passed <- NA
      result$Runs_p <- NA
      result$Bias_passed <- NA
    }
    
    # 3. Parameter estimates (terminal year)
    result$BBmsy_terminal <- NA
    result$FFmsy_terminal <- NA
    
    if (!is.null(jb$kobe) && is.data.frame(jb$kobe) && nrow(jb$kobe) > 0) {
      terminal_kobe <- jb$kobe[nrow(jb$kobe), ]
      if ("stock" %in% names(terminal_kobe)) {
        result$BBmsy_terminal <- terminal_kobe$stock[1]
      }
      if ("harvest" %in% names(terminal_kobe)) {
        result$FFmsy_terminal <- terminal_kobe$harvest[1]
      }
    }
    
    return(result)
    
  }, error = function(e) {
    # Return minimal result with NAs on error
    warning("Error in jabbaDiagnostics for ", id, " - ", scenario, ": ", e$message)
    result <- data.frame(
      Stock = if (is.null(id)) NA_character_ else id,
      Scenario = if (is.null(scenario)) NA_character_ else scenario,
      Geweke_passed = NA,
      Geweke_mean = NA,
      Heidel_passed = NA,
      Heidel_mean = NA,
      n_residuals = NA,
      mean_residual = NA,
      sd_residual = NA,
      Normality_passed = NA,
      Normality_p = NA,
      Autocorr_passed = NA,
      Autocorr_p = NA,
      Runs_passed = NA,
      Runs_p = NA,
      Bias_passed = NA,
      BBmsy_terminal = NA,
      FFmsy_terminal = NA,
      stringsAsFactors = FALSE
    )
    return(result)
  })
}

#' Summarize Diagnostics Across Multiple Fits
#'
#' @description
#' Creates a summary table of diagnostic results across multiple JABBA fits,
#' grouped by scenario. Useful for reporting diagnostic pass rates.
#'
#' @param diagnostics A data.frame from \code{jabbaDiagnostics()} containing
#'   diagnostic results for multiple fits
#' @param scenarios Character vector. Scenarios to include in summary.
#'   If NULL, includes all scenarios. Default: NULL
#' @param scenario_names Named character vector. Optional mapping from scenario
#'   names to display names (e.g., c("EBiomass" = "JABBA-ICES")). Default: NULL
#'
#' @return A data.frame with summary statistics by scenario:
#'   \item{Scenario_name}{Display name of scenario}
#'   \item{n_stocks}{Number of stocks}
#'   \item{n_converged_geweke}{Number passing Geweke test}
#'   \item{n_converged_heidel}{Number passing Heidel test}
#'   \item{n_normal}{Number passing normality test}
#'   \item{n_normal_total}{Total with normality test results}
#'   \item{n_no_autocorr}{Number passing autocorrelation test}
#'   \item{n_no_autocorr_total}{Total with autocorrelation test results}
#'   \item{n_runs_passed}{Number passing runs test}
#'   \item{n_runs_total}{Total with runs test results}
#'   \item{n_bias_passed}{Number passing bias test}
#'   \item{n_bias_total}{Total with bias test results}
#'
#' @examples
#' \dontrun{
#' # Run diagnostics on all fits
#' diags <- jabbaDiagnostics(fits)
#'
#' # Summarize by scenario
#' summary <- summarizeDiagnostics(diags, 
#'                                 scenarios = c("EBiomass", "EBiomass 1903"),
#'                                 scenario_names = c("EBiomass" = "JABBA-ICES",
#'                                                   "EBiomass 1903" = "JABBA-Historical"))
#' }
#'
#' @export
summarizeDiagnostics <- function(diagnostics, scenarios = NULL, scenario_names = NULL) {
  
  if (is.null(diagnostics) || nrow(diagnostics) == 0) {
    stop("diagnostics data.frame is empty or NULL")
  }
  
  # Filter scenarios if specified
  if (!is.null(scenarios)) {
    diagnostics <- diagnostics[diagnostics$Scenario %in% scenarios, ]
  }
  
  if (nrow(diagnostics) == 0) {
    stop("No data remaining after filtering scenarios")
  }
  
  # Summarize by scenario
  summary_by_scenario <- diagnostics %>%
    dplyr::group_by(.data$Scenario) %>%
    dplyr::summarize(
      n_stocks = dplyr::n(),
      n_converged_geweke = sum(.data$Geweke_passed, na.rm = TRUE),
      n_converged_heidel = sum(.data$Heidel_passed, na.rm = TRUE),
      n_normal = sum(.data$Normality_passed, na.rm = TRUE),
      n_normal_total = sum(!is.na(.data$Normality_passed)),
      n_no_autocorr = sum(.data$Autocorr_passed, na.rm = TRUE),
      n_no_autocorr_total = sum(!is.na(.data$Autocorr_passed)),
      n_runs_passed = sum(.data$Runs_passed, na.rm = TRUE),
      n_runs_total = sum(!is.na(.data$Runs_passed)),
      n_bias_passed = sum(.data$Bias_passed, na.rm = TRUE),
      n_bias_total = sum(!is.na(.data$Bias_passed)),
      .groups = "drop"
    )
  
  # Apply scenario name mapping if provided
  if (!is.null(scenario_names)) {
    summary_by_scenario$Scenario_name <- scenario_names[summary_by_scenario$Scenario]
    # Use original name if mapping not found
    missing <- is.na(summary_by_scenario$Scenario_name)
    summary_by_scenario$Scenario_name[missing] <- summary_by_scenario$Scenario[missing]
  } else {
    summary_by_scenario$Scenario_name <- summary_by_scenario$Scenario
  }
  
  return(summary_by_scenario)
}
