#' Extract JABBA Time Series as FLQuants
#'
#' @description
#' Converts selected slices of a JABBA \code{fit$timeseries} array into an
#' \code{FLQuants} object for FLR plotting and process-error diagnostics.
#'
#' @param fit A JABBA fit with a 3-d \code{timeseries} array
#'   (year x statistic x quantity). Also accepts wrappers with
#'   \code{fit$fit$timeseries}.
#' @param quant Character. Second-dimension statistic. Default \code{"mu"}.
#' @param vars Named character vector mapping output names to JABBA quantities.
#'   Default maps \code{B}, \code{F}, \code{BBmsy}, \code{FFmsy}, \code{procB},
#'   \code{SPt} to \code{stock}, \code{harvest}, \code{bbmsy}, \code{ffmsy},
#'   \code{pe}, \code{sprod}. Missing names are skipped with a warning.
#'
#' @return An \code{FLQuants} of annual series.
#'
#' @details
#' JABBA stores process error as \code{procB}. This exposes that series
#' directly; see \code{\link{jabbaPE}} for production-function residuals.
#'
#' @examples
#' \dontrun{
#' fit <- runJABBA(stk, method = "ices", quick = TRUE)
#' ts  <- jabbaTs(fit)
#' plotPe(ts)
#' }
#'
#' @seealso \code{\link{jabbaTsCI}}, \code{\link{plotJabbaTs}},
#'   \code{\link{plotPe}}, \code{\link{jabbaPE}}, \code{\link{rod}}
#' @export
#' @importFrom FLCore FLQuants
jabbaTs <- function(fit,
                    quant = "mu",
                    vars = c(stock = "B",
                             harvest = "F",
                             bbmsy = "BBmsy",
                             ffmsy = "FFmsy",
                             pe = "procB",
                             sprod = "SPt")) {

  u <- .jabbaTsArray(fit)
  if (!quant %in% u$stats)
    stop("Statistic '", quant, "' not found in fit$timeseries. Available: ",
         paste(u$stats, collapse = ", "))

  out <- list()
  for (nm in names(vars)) {
    jnm <- vars[[nm]]
    if (!jnm %in% u$quants) {
      warning("Quantity '", jnm, "' (as '", nm,
              "') not found in fit$timeseries; skipping")
      next
    }
    val <- as.numeric(u$ts[, quant, jnm])
    if (length(val) != length(u$years))
      stop("timeseries '", jnm, "' length does not match years")
    out[[nm]] <- FLCore::FLQuant(val, dimnames = list(year = u$years))
  }

  if (!length(out))
    stop("No requested timeseries quantities were found in fit$timeseries")

  do.call(FLCore::FLQuants, out)
}

.jabbaTsArray <- function(fit) {
  if (is.null(fit))
    stop("'fit' cannot be NULL")
  if (is.null(fit$timeseries) && !is.null(fit$fit$timeseries))
    fit <- fit$fit
  ts <- fit$timeseries
  if (is.null(ts))
    stop("JABBA fit object must have a 'timeseries' component")
  if (length(dim(ts)) != 3)
    stop("'fit$timeseries' must be a 3-d array (year x statistic x quantity)")
  dn <- dimnames(ts)
  years <- as.numeric(dn[[1]])
  if (!length(years) || any(!is.finite(years)))
    stop("fit$timeseries year dimnames must be numeric")
  list(ts = ts, years = years, stats = dn[[2]], quants = dn[[3]])
}

#' Median and 95% CI from a JABBA fit
#'
#' Long data frame of \code{mu}, \code{lci} and \code{uci} from
#' \code{fit$timeseries}, for ggplot ribbons.
#'
#' @param fit A JABBA fit (or list with \code{fit$timeseries}).
#' @param vars Named character vector as in \code{\link{jabbaTs}}.
#'
#' @return A data frame with \code{year}, \code{qname}, \code{mu}, \code{lci},
#'   \code{uci}.
#'
#' @examples
#' \dontrun{
#' ci <- jabbaTsCI(fit)
#' ggplot(ci, aes(year)) +
#'   geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
#'   geom_line(aes(y = mu)) +
#'   facet_wrap(~qname, scales = "free_y")
#' }
#'
#' @seealso \code{\link{jabbaTs}}, \code{\link{plotJabbaTs}}
#' @export
jabbaTsCI <- function(fit,
                     vars = c(stock = "B",
                              harvest = "F",
                              bbmsy = "BBmsy",
                              ffmsy = "FFmsy",
                              pe = "procB",
                              sprod = "SPt")) {
  u <- .jabbaTsArray(fit)
  need <- c("mu", "lci", "uci")
  if (!all(need %in% u$stats))
    stop("fit$timeseries must have statistics mu, lci, uci")

  rows <- list()
  for (nm in names(vars)) {
    jnm <- vars[[nm]]
    if (!jnm %in% u$quants) {
      warning("Quantity '", jnm, "' (as '", nm,
              "') not found in fit$timeseries; skipping")
      next
    }
    rows[[nm]] <- data.frame(
      year = u$years,
      qname = nm,
      mu = as.numeric(u$ts[, "mu", jnm]),
      lci = as.numeric(u$ts[, "lci", jnm]),
      uci = as.numeric(u$ts[, "uci", jnm]),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows))
    stop("No requested timeseries quantities were found in fit$timeseries")
  do.call(rbind, rows)
}

#' Plot JABBA trajectories with 95% CI ribbons
#'
#' @param x A JABBA fit, or a named list of fits to overlay.
#' @param vars Named character vector as in \code{\link{jabbaTs}}.
#' @param qname Optional subset of \code{jabbaTsCI()} names (e.g.
#'   \code{"bbmsy"}).
#' @param om Optional OM overlay: an \code{FLQuant}, a named list of
#'   \code{FLQuant}s, or \code{FLQuants}. Typically OM \(B/B_{\mathrm{MSY}}\)
#'   (SSB and/or exploitable biomass). Drawn on the \code{om.qname} panel.
#' @param om.qname Panel that receives \code{om}. Default \code{"bbmsy"}.
#' @param ... Unused.
#'
#' @return A ggplot.
#'
#' @examples
#' \dontrun{
#' plotJabbaTs(fit)
#' plotJabbaTs(list(JABBA = fit, `JABBA-PE` = fitPE), qname = "bbmsy")
#' }
#'
#' @seealso \code{\link{jabbaTsCI}}, \code{\link{jabbaTs}}
#' @export
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line facet_wrap labs
#'   theme_bw
plotJabbaTs <- function(x,
                        vars = c(stock = "B",
                                 harvest = "F",
                                 bbmsy = "BBmsy",
                                 ffmsy = "FFmsy",
                                 pe = "procB",
                                 sprod = "SPt"),
                        qname = NULL,
                        om = NULL, om.qname = "bbmsy", ...) {
  if (!is.null(qname)) {
    keep <- names(vars) %in% qname
    if (!any(keep))
      stop("qname not found in vars: ", paste(qname, collapse = ", "))
    vars <- vars[keep]
  }
  is.fit <- is.list(x) && (!is.null(x$timeseries) || !is.null(x$fit$timeseries))
  if (is.fit) {
    d <- cbind(jabbaTsCI(x, vars = vars), run = "JABBA", stringsAsFactors = FALSE)
  } else {
    if (is.null(names(x)) || any(!nzchar(names(x))))
      stop("A list of fits must be named, e.g. list(JABBA = fit, `JABBA-PE` = fitPE)")
    d <- do.call(rbind, lapply(names(x), function(nm) {
      cbind(jabbaTsCI(x[[nm]], vars = vars), run = nm, stringsAsFactors = FALSE)
    }))
  }
  if (!is.null(qname))
    d <- d[d$qname %in% qname, , drop = FALSE]
  if (!nrow(d))
    stop("No series to plot")

  nq <- length(unique(d$qname))
  p <- ggplot2::ggplot(d, ggplot2::aes(year, fill = run, colour = run)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci),
                          alpha = 0.25, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = mu)) +
    ggplot2::labs(x = "Year", y = NULL, colour = NULL, fill = NULL) +
    ggplot2::theme_bw()
  if (nq > 1L)
    p <- p + ggplot2::facet_wrap(~qname, scales = "free_y", ncol = 1)

  if (!is.null(om)) {
    omd <- .jabbaOmDf(om, om.qname = om.qname, qname = qname)
    if (!is.null(omd) && nrow(omd))
      p <- p + ggplot2::geom_line(
        data = omd,
        ggplot2::aes(year, mu, colour = run),
        inherit.aes = FALSE,
        linetype = 2
      )
  }
  p
}

.jabbaOmDf <- function(om, om.qname = "bbmsy", qname = NULL) {
  if (is.null(om))
    return(NULL)
  if (methods::is(om, "FLQuant"))
    om <- list(OM = om)
  else if (!(methods::is(om, "FLQuants") || is.list(om)))
    stop("'om' must be an FLQuant, FLQuants, or named list of FLQuants")
  nms <- names(om)
  if (is.null(nms) || any(!nzchar(nms)))
    stop("A list of OM series must be named, e.g. FLQuants(`OM SSB` = ..., `OM EB` = ...)")
  rows <- lapply(nms, function(nm) {
    x <- om[[nm]]
    if (!methods::is(x, "FLQuant"))
      stop("om['", nm, "'] must be an FLQuant")
    data.frame(
      year = as.numeric(dimnames(x)$year),
      mu = as.numeric(x),
      qname = om.qname,
      run = nm,
      stringsAsFactors = FALSE
    )
  })
  omd <- do.call(rbind, rows)
  if (!is.null(qname))
    omd <- omd[omd$qname %in% qname, , drop = FALSE]
  omd
}
