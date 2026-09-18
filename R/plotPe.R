#' Four-Panel Process-Error Diagnostics
#'
#' @description
#' Process-error residual diagnostics: (1) residuals vs time with optional
#' STARS regime polygons, (2) residuals vs relative biomass, (3) ACF,
#' (4) residual histogram with normal density overlay.
#'
#' @param object \code{FLQuants} with at least \code{pe} (preferably also
#'   \code{bbmsy}); or an \code{FLQuant} of residuals; or a JABBA fit
#'   (passed through \code{\link{jabbaTs}}).
#' @param pe Name of the process-error series. Default \code{"pe"}.
#' @param status Name of the relative biomass series for panel 2.
#'   Default \code{"bbmsy"}; falls back to \code{"stock"}.
#' @param lag.max Maximum ACF lag (default 10).
#' @param bins Histogram bins (default 40).
#' @param regimes Overlay STARS polygons via \code{\link{rod}} (default TRUE).
#' @param n,sig Arguments forwarded to \code{\link{rod}}.
#' @param combine If TRUE and \pkg{patchwork} is available, return a 2x2
#'   combined plot; otherwise a named list of ggplots.
#' @param ... Unused.
#'
#' @return A patchwork object, or a named list with elements
#'   \code{time}, \code{status}, \code{acf}, \code{hist}.
#'
#' @examples
#' \dontrun{
#' plotPe(jabbaTs(fit))
#' }
#'
#' @seealso \code{\link{jabbaTs}}, \code{\link{rod}}, \code{\link{jabbaPE}}
#' @export
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_point geom_smooth
#'   geom_col geom_histogram geom_polygon geom_line labs theme_minimal
#'   stat_function after_stat
plotPe <- function(object,
                   pe = "pe",
                   status = "bbmsy",
                   lag.max = 10,
                   bins = 40,
                   regimes = TRUE,
                   n = 10,
                   sig = 1.68,
                   combine = TRUE,
                   ...) {

  if (inherits(object, "FLQuant")) {
    object <- FLCore::FLQuants(pe = object)
    pe <- "pe"
  } else if (is.list(object) && (!is.null(object$timeseries) ||
                                  !is.null(object$fit$timeseries))) {
    object <- jabbaTs(object)
  }

  if (!inherits(object, "FLQuants"))
    stop("'object' must be FLQuants, FLQuant, or a JABBA fit")
  if (!pe %in% names(object))
    stop("Process-error series '", pe, "' not found in object")

  peQ <- object[[pe]]
  peVals <- as.numeric(c(peQ))
  peDF <- as.data.frame(peQ)
  peDF$year <- as.numeric(as.character(peDF$year))

  ## STARS regime bands first (background)
  p1 <- ggplot2::ggplot(peDF, ggplot2::aes(x = year, y = data)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = 2)

  if (isTRUE(regimes)) {
    stars <- try(
      rodFn(data = as.numeric(peDF$data),
            year = as.numeric(peDF$year),
            n = n, sig = sig, plot = FALSE),
      silent = TRUE
    )
    if (inherits(stars, "try-error")) {
      warning("STARS regime detection failed: ",
              as.character(stars))
    } else if (is.data.frame(stars) && nrow(stars)) {
      rects <- unique(stars[, c("regime", "minyear", "maxyear", "mn", "sd")])
      rects$ymin <- rects$mn - rects$sd
      rects$ymax <- rects$mn + rects$sd
      p1 <- p1 +
        ggplot2::geom_rect(
          data = rects,
          ggplot2::aes(xmin = minyear, xmax = maxyear,
                       ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill = "steelblue", alpha = 0.25, colour = NA
        )
    }
  }

  p1 <- p1 +
    ggplot2::geom_line(colour = "grey30") +
    ggplot2::geom_point(alpha = 0.7, colour = "steelblue") +
    ggplot2::labs(x = "Year", y = "Process residual") +
    ggplot2::theme_minimal()

  xnm <- if (status %in% names(object)) {
    status
  } else if ("stock" %in% names(object)) {
    "stock"
  } else {
    NA_character_
  }

  if (is.na(xnm)) {
    p2 <- ggplot2::ggplot() +
      ggplot2::labs(title = "No status / stock series for residual plot") +
      ggplot2::theme_minimal()
  } else {
    ## Bind columns explicitly (do not use stats::model.frame).
    ## Use status/resid names so formula y ~ x is unambiguous for geom_smooth.
    pe_df <- as.data.frame(object[[pe]])
    st_df <- as.data.frame(object[[xnm]])
    mf <- data.frame(
      year  = as.numeric(pe_df$year),
      resid = as.numeric(pe_df$data),
      status = as.numeric(st_df$data)
    )
    xlab <- if (identical(xnm, "bbmsy")) expression(B/B[MSY]) else xnm
    p2 <- ggplot2::ggplot(mf, ggplot2::aes(status, resid)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
      ggplot2::geom_point(alpha = 0.6, colour = "steelblue") +
      ggplot2::geom_smooth(se = FALSE, colour = "darkorange", method = "loess") +
      ggplot2::labs(x = xlab, y = "Process residual") +
      ggplot2::theme_minimal()
    if (identical(xnm, "bbmsy"))
      p2 <- p2 + ggplot2::geom_vline(xintercept = 1, colour = "red",
                                       linetype = 2)
  }

  ac <- stats::acf(peVals, lag.max = lag.max, plot = FALSE,
                   na.action = stats::na.pass)
  acfDF <- data.frame(lag = as.numeric(ac$lag), acf = as.numeric(ac$acf))
  crit <- stats::qnorm(0.975) / sqrt(sum(is.finite(peVals)))
  p3 <- ggplot2::ggplot(acfDF, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_hline(yintercept = c(-crit, crit), linetype = "dashed",
                        colour = "grey60") +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(x = "Lag", y = "ACF") +
    ggplot2::theme_minimal()

  mu <- mean(peVals, na.rm = TRUE)
  sdv <- stats::sd(peVals, na.rm = TRUE)
  p4 <- ggplot2::ggplot(peDF, ggplot2::aes(x = data)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = bins, fill = "steelblue",
                            colour = "white", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = 0, colour = "red", linetype = 2) +
    ggplot2::stat_function(
      fun = stats::dnorm,
      args = list(mean = mu, sd = sdv),
      colour = "red", linewidth = 0.8
    ) +
    ggplot2::labs(x = "Residual", y = "Density") +
    ggplot2::theme_minimal()

  
  p1=p1+labs(subtitle="Time Series")
  p2=p2+labs(subtitle="Depletion")
  p3=p3+labs(subtitle="Lags")
  p4=p4+labs(subtitle="Mean and Variance")
  plots <- list(time = p1, status = p2, acf = p3, hist = p4)

  if (isTRUE(combine) && requireNamespace("patchwork", quietly = TRUE))
    return((plots$time + plots$status) / (plots$acf + plots$hist))

  if (isTRUE(combine))
    warning("Package 'patchwork' not installed; returning list of ggplots")

  plots
}
