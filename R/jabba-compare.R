#' JABBA depletion-dependent process-error scale
#'
#' Relative process SD used by JABBA when \code{proc.dyn = TRUE}.
#' Normalised so \code{relsig(1) = 1} at \eqn{B/B_{\mathrm{MSY}}}.
#' Defaults are the ICES PE-SD coefficients hardcoded in JABBA
#' (\code{a = 0.128}, \code{c = 0.223}, \code{d = 1.326}). Override from
#' \code{\link{runJABBA}} with \code{pe = c(a = 0.2, c = 0.35, d = 2)}.
#'
#' \deqn{\mathrm{relsig}(x)=\frac{a+(1-a)/(1+(x/c)^{d})}{a+(1-a)/(1+(1/c)^{d})}}
#'
#' @param x Relative biomass \eqn{B/B_{\mathrm{MSY}}}.
#' @param a Floor as \eqn{x \to \infty}. Default \code{0.128}.
#' @param c Inflection (where variance rises). Default \code{0.223}.
#' @param d Steepness of that rise. Default \code{1.326}.
#'
#' @return Numeric vector the same length as \code{x}.
#'
#' @examples
#' relsig(1)
#' relsig(c(0.2, 1, 2))
#'
#' @seealso \code{\link{compareJabbaPE}}, \code{\link{runJABBA}}
#' @export
relsig <- function(x, a = 0.128, c = 0.223, d = 1.326) {
  raw <- function(z) a + (1 - a) / (1 + (z / c)^d)
  raw(x) / raw(1)
}

.unwrapJabbaFit <- function(x) {
  if (is.null(x))
    return(NULL)
  if (!is.null(x$timeseries))
    return(x)
  if (!is.null(x$fit$timeseries))
    return(x$fit)
  NULL
}

.jabbaFitOk <- function(x) {
  !is.null(.unwrapJabbaFit(x))
}

.jabbaParMedian <- function(fit, par) {
  f <- .unwrapJabbaFit(fit)
  if (is.null(f) || is.null(f$pars))
    return(NA_real_)
  p <- f$pars
  rn <- rownames(p)
  col <- if ("Median" %in% colnames(p)) "Median" else colnames(p)[1]
  if (!is.null(rn) && par %in% rn)
    return(as.numeric(p[par, col]))
  NA_real_
}

.jabbaParsTable <- function(fit, run) {
  f <- .unwrapJabbaFit(fit)
  if (is.null(f) || is.null(f$pars))
    return(NULL)
  p <- as.data.frame(f$pars, stringsAsFactors = FALSE)
  p$parameter <- rownames(f$pars)
  rownames(p) <- NULL
  p$run <- run
  p[, c("run", "parameter", setdiff(names(p), c("run", "parameter")))]
}

.jabbaStatusRow <- function(fit, run) {
  f <- .unwrapJabbaFit(fit)
  ts <- try(jabbaTs(fit), silent = TRUE)
  last <- function(q) {
    v <- as.numeric(c(q))
    v[length(v)]
  }
  bb <- ff <- NA_real_
  if (!inherits(ts, "try-error")) {
    if ("bbmsy" %in% names(ts))
      bb <- last(ts$bbmsy)
    if ("ffmsy" %in% names(ts))
      ff <- last(ts$ffmsy)
  }
  data.frame(
    run = run,
    BBmsy = bb,
    FFmsy = ff,
    sigma.proc = .jabbaParMedian(fit, "sigma.proc"),
    r = .jabbaParMedian(fit, "r"),
    K = .jabbaParMedian(fit, "K"),
    converged = .jabbaFitOk(fit),
    stringsAsFactors = FALSE
  )
}

.compareJabbaCall <- function(object, eq, proc.dyn, scenario, dots) {
  dots$output <- NULL
  args <- c(list(object = object, proc.dyn = proc.dyn, scenario = scenario,
                 output = "jabba"), dots)
  if (methods::is(object, "FLStock") && !is.null(eq))
    args$eq <- eq
  do.call(runJABBA, args)
}

#' Compare JABBA with and without depletion-dependent process error
#'
#' Fits the same catch and biomass index twice: constant process SD
#' (\code{proc.dyn = FALSE}) and JABBA's depletion-dependent scaling
#' (\code{proc.dyn = TRUE}), where
#' \eqn{\sigma_t = \sigma \times \mathrm{relsig}(B_t/B_{\mathrm{MSY}})}.
#'
#' Requires a JABBA install that exposes \code{proc.dyn} in
#' \code{build_jabba()} (current GitHub \code{jabbamodel/JABBA}).
#' The \code{a}, \code{c}, \code{d} coefficients are those hardcoded in
#' JABBA, not estimated here; see \code{\link{relsig}}.
#'
#' @param object An \code{FLStock} or catch \code{data.frame}, as for
#'   \code{\link{runJABBA}}.
#' @param eq Optional \code{FLBRP} used for production-function priors.
#' @param ... Passed to both \code{runJABBA()} calls (priors, index,
#'   \code{quick}, \code{nc}, \ldots). \code{proc.dyn} in \code{...} is
#'   ignored so the two arms stay distinct.
#'
#' @return An object of class \code{compareJabbaPE}:
#' \itemize{
#'   \item \code{fits}: named list \code{constant}, \code{proc.dyn}
#'   \item \code{ts}: \code{FLQuants} from \code{\link{jabbaTs}} for each arm
#'   \item \code{pars}: combined JABBA parameter table
#'   \item \code{status}: terminal \eqn{B/B_{\mathrm{MSY}}}, \eqn{F/F_{\mathrm{MSY}}},
#'     and \code{sigma.proc}
#'   \item \code{relsig}: the scaling curve used by JABBA
#' }
#'
#' @examples
#' \dontrun{
#' data(ple4, package = "FLCore")
#' data(ple4brp, package = "FLBRP")
#' cmp <- compareJabbaPE(ple4, eq = ple4brp, quick = TRUE)
#' cmp$status
#' plot(cmp)
#' plotPe(cmp$ts$constant)
#' plotPe(cmp$ts$proc.dyn)
#' }
#'
#' @seealso \code{\link{runJABBA}}, \code{\link{relsig}}, \code{\link{plotPe}},
#'   \code{\link{jabbaTs}}
#' @export
#' @importFrom methods setOldClass
compareJabbaPE <- function(object, eq = NULL, ...) {
  if (!.jabbaHasProcDyn())
    stop("Installed JABBA has no proc.dyn argument. Update with ",
         "remotes::install_github(\"jabbamodel/JABBA\")")

  dots <- list(...)
  dots$proc.dyn <- NULL
  base <- dots$scenario
  dots$scenario <- NULL
  if (is.null(base) || !nzchar(as.character(base)[1])) {
    scen0 <- "constant"
    scen1 <- "proc.dyn"
  } else {
    scen0 <- paste0(base, "_constant")
    scen1 <- paste0(base, "_proc.dyn")
  }

  fit0 <- .compareJabbaCall(object, eq, proc.dyn = FALSE, scenario = scen0, dots)
  fit1 <- .compareJabbaCall(object, eq, proc.dyn = TRUE,  scenario = scen1, dots)

  if (!.jabbaFitOk(fit0))
    warning("Constant-PE JABBA fit did not return a timeseries")
  if (!.jabbaFitOk(fit1))
    warning("Depletion-dependent (proc.dyn) JABBA fit did not return a timeseries")

  ts <- list()
  if (.jabbaFitOk(fit0))
    ts$constant <- jabbaTs(fit0)
  if (.jabbaFitOk(fit1))
    ts$proc.dyn <- jabbaTs(fit1)

  pars <- rbind(
    .jabbaParsTable(fit0, "constant"),
    .jabbaParsTable(fit1, "proc.dyn")
  )
  status <- rbind(
    .jabbaStatusRow(fit0, "constant"),
    .jabbaStatusRow(fit1, "proc.dyn")
  )

  x <- seq(0.02, 3, length.out = 200)
  rs <- do.call(.jabbaRelsigPars, dots)
  curve <- data.frame(bbmsy = x, relsig = relsig(x, a = rs$a, c = rs$c, d = rs$d))

  out <- list(
    fits = list(constant = fit0, proc.dyn = fit1),
    ts = ts,
    pars = pars,
    status = status,
    relsig = curve,
    relsig.par = c(a = rs$a, c = rs$c, d = rs$d)
  )
  class(out) <- c("compareJabbaPE", "list")
  out
}

methods::setOldClass("compareJabbaPE")

#' @export
#' @rdname compareJabbaPE
print.compareJabbaPE <- function(x, ...) {
  cat("JABBA process-error comparison (constant vs proc.dyn)\n")
  ok0 <- .jabbaFitOk(x$fits$constant)
  ok1 <- .jabbaFitOk(x$fits$proc.dyn)
  cat("  constant : ", if (ok0) "fitted" else "failed", "\n", sep = "")
  cat("  proc.dyn : ", if (ok1) "fitted" else "failed", "\n", sep = "")
  if (!is.null(x$status)) {
    cat("\nTerminal status:\n")
    print(x$status, row.names = FALSE)
  }
  invisible(x)
}

#' Plot constant vs depletion-dependent JABBA fits
#'
#' Overlay \eqn{B/B_{\mathrm{MSY}}}, \eqn{F/F_{\mathrm{MSY}}} and process
#' residuals, plus the \code{\link{relsig}} curve used by \code{proc.dyn}.
#'
#' @param x A \code{compareJabbaPE} object.
#' @param ... Unused.
#'
#' @return A patchwork object if \pkg{patchwork} is available, otherwise a
#'   named list of ggplots.
#' @export
#' @rdname compareJabbaPE
plot.compareJabbaPE <- function(x, ...) {
  rows <- list()
  for (run in names(x$ts)) {
    ts <- x$ts[[run]]
    df <- as.data.frame(ts)
    df$year <- as.numeric(as.character(df$year))
    df$run <- run
    rows[[run]] <- df
  }
  if (!length(rows))
    stop("No fitted timeseries to plot")
  df <- do.call(rbind, rows)
  keep <- intersect(c("bbmsy", "ffmsy", "pe"), unique(df$qname))
  df <- df[df$qname %in% keep, , drop = FALSE]
  lab <- c(bbmsy = "B/Bmsy", ffmsy = "F/Fmsy", pe = "Process residual")
  df$panel <- lab[as.character(df$qname)]
  df$panel[is.na(df$panel)] <- as.character(df$qname[is.na(df$panel)])

  p.ts <- ggplot2::ggplot(df, ggplot2::aes(year, data, colour = run)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linetype = 2) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "Year", y = NULL, colour = NULL) +
    ggplot2::theme_minimal()

  p.sig <- ggplot2::ggplot(x$relsig, ggplot2::aes(bbmsy, relsig)) +
    ggplot2::geom_hline(yintercept = 1, colour = "grey70", linetype = 2) +
    ggplot2::geom_vline(xintercept = 1, colour = "red", linetype = 2) +
    ggplot2::geom_line(colour = "grey20") +
    ggplot2::labs(x = expression(B/B[MSY]),
                  y = expression(sigma/sigma[MSY]),
                  title = "proc.dyn scale") +
    ggplot2::theme_minimal()

  if ("proc.dyn" %in% names(x$ts) && "bbmsy" %in% names(x$ts$proc.dyn)) {
    rp <- x$relsig.par
    if (is.null(rp))
      rp <- c(a = 0.128, c = 0.223, d = 1.326)
    bb <- as.data.frame(x$ts$proc.dyn$bbmsy)
    p.sig <- p.sig +
      ggplot2::geom_point(
        data = data.frame(
          bbmsy = as.numeric(bb$data),
          relsig = relsig(as.numeric(bb$data),
                           a = unname(rp[["a"]]),
                           c = unname(rp[["c"]]),
                           d = unname(rp[["d"]]))),
        colour = "steelblue", alpha = 0.5
      )
  }

  if (requireNamespace("patchwork", quietly = TRUE))
    return(p.ts | p.sig)

  warning("Package 'patchwork' not installed; returning list of ggplots")
  list(timeseries = p.ts, relsig = p.sig)
}

#' @rdname compareJabbaPE
#' @export
setMethod("plot", signature(x = "compareJabbaPE"),
          function(x, y, ...) plot.compareJabbaPE(x, ...))

#' @rdname compareJabbaPE
#' @export
setMethod("print", signature(x = "compareJabbaPE"),
          function(x, ...) print.compareJabbaPE(x, ...))

#' @rdname compareJabbaPE
#' @export
setMethod("show", signature(object = "compareJabbaPE"),
          function(object) print.compareJabbaPE(object))

