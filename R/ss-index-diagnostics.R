#' Stock Synthesis index standardization and diagnostics
#'
#' Standardize index series by group, fit smooth trends, and derive residual
#' diagnostics similar to the CPUE workflow used in the WKBSEABASS scripts.
#'
#' Prepare SS index data for diagnostics
#'
#' Convert `cpueSS`-style list output into a standard data.frame with
#' `name`, `year`, and `data` columns.
#'
#' @param object A list-like object such as output from `cpueSS`.
#'
#' @return For `ssIndexPrep`, a data.frame with `name`, `year`, and `data`.
#' @export
setMethod("ssIndexPrep", signature(object = "list"), function(object) {
  dat <- plyr::ldply(object, as.data.frame, drop = TRUE)
  if (!all(c(".id", "year", "data") %in% names(dat))) {
    stop("Expected columns '.id', 'year', and 'data' in list-expanded cpue data.")
  }
  dat <- dat[!is.na(dat$data) & !is.na(dat$year), c(".id", "year", "data"), drop = FALSE]
  names(dat) <- c("name", "year", "data")
  dat
})

#' Prepare SS index data from a data.frame
#'
#' @param object A data.frame containing index observations.
#' @param name Name of the index/series column.
#' @param year Name of the year column.
#' @param data Name of the observed index value column.
#'
#' @return For `ssIndexPrep`, a data.frame with `name`, `year`, and `data`.
#' @export
setMethod("ssIndexPrep", signature(object = "data.frame"), function(
  object, name = "name", year = "year", data = "data"
) {
  stopifnot(all(c(name, year, data) %in% names(object)))
  out <- object[, c(name, year, data), drop = FALSE]
  names(out) <- c("name", "year", "data")
  out <- out[!is.na(out$data) & !is.na(out$year), , drop = FALSE]
  out
})

#' Standardize SS index data and fit smooth trends
#'
#' @param object A data.frame containing at least year, value, and group fields.
#' @param year Name of the year column.
#' @param value Name of the observed index value column.
#' @param group Name of the grouping column (index/series name).
#' @param span Smoothing span passed to `stats::lowess`.
#' @param center_zero Logical; if `TRUE`, standardize each group to mean 0 and sd 1.
#'
#' @return For `ssIndexScale`, a data.frame with columns `name`, `year`, `data`,
#'   `hat`, and `residual`.
#' @export
setMethod("ssIndexScale", signature(object = "data.frame"), function(
  object, year = "year", value = "data", group = "name", span = 2 / 3, center_zero = TRUE
) {
  stopifnot(all(c(year, value, group) %in% names(object)))

  dat <- data.frame(
    name = object[[group]],
    year = as.numeric(object[[year]]),
    data = as.numeric(object[[value]])
  )
  dat <- dat[!is.na(dat$year) & !is.na(dat$data) & !is.na(dat$name), , drop = FALSE]

  split_dat <- split(dat, dat$name)
  scaled <- lapply(split_dat, function(df) {
    if (center_zero) {
      s <- stats::sd(df$data, na.rm = TRUE)
      if (is.na(s) || s == 0) s <- 1
      df$data <- (df$data - mean(df$data, na.rm = TRUE)) / s
    }
    sm <- stats::lowess(df$year, df$data, f = span)
    df$hat <- stats::approx(sm$x, sm$y, xout = df$year, ties = mean)$y
    df$residual <- df$data - df$hat
    df
  })

  out <- do.call(rbind, scaled)
  rownames(out) <- NULL
  out[order(out$name, out$year), c("name", "year", "data", "hat", "residual")]
})

#' Runs test for standardized SS index residuals
#'
#' @param object A data.frame as returned by `ssIndexScale`.
#' @param alpha Significance level for the two-sided runs test.
#'
#' @return For `ssIndexRuns`, a data.frame with runs-test diagnostics by index.
#' @export
setMethod("ssIndexRuns", signature(object = "data.frame"), function(object, alpha = 0.05) {
  stopifnot(all(c("name", "residual") %in% names(object)))

  run_test <- function(x) {
    x <- x[!is.na(x) & x != 0]
    n <- length(x)
    if (n < 3) {
      return(data.frame(runs = NA_integer_, z = NA_real_, p_value = NA_real_, pass = NA))
    }
    s <- ifelse(x > 0, 1L, -1L)
    r <- length(rle(s)$lengths)
    n1 <- sum(s > 0)
    n2 <- sum(s < 0)
    mu <- (2 * n1 * n2) / (n1 + n2) + 1
    vr <- (2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / (((n1 + n2)^2) * (n1 + n2 - 1))
    z <- if (vr > 0) (r - mu) / sqrt(vr) else NA_real_
    p <- if (!is.na(z)) 2 * (1 - stats::pnorm(abs(z))) else NA_real_
    data.frame(runs = r, z = z, p_value = p, pass = ifelse(is.na(p), NA, p > alpha))
  }

  by_name <- split(object$residual, object$name)
  out <- lapply(names(by_name), function(nm) cbind(name = nm, run_test(by_name[[nm]])))
  do.call(rbind, out)
})

#' Build combined SS index diagnostics
#'
#' @param object A data.frame containing index observations.
#' @param ... Additional arguments passed to `ssIndexScale` and `ssIndexRuns`.
#'
#' @return For `ssIndexDiagnostics`, a list with `scaled` and `runs`.
#' @export
setMethod("ssIndexDiagnostics", signature(object = "data.frame"), function(object, ...) {
  scaled <- ssIndexScale(object, ...)
  runs <- ssIndexRuns(scaled)
  list(scaled = scaled, runs = runs)
})
