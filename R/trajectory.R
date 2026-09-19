# trajectory.R - relative biomass and fishing mortality trajectories
# FLRebuild
#
# S4 generics and methods, moved from 09_plot_ssf_from_report.R.
#
# Companion to mvln.R. Both treat the SS3 point estimate as the median
# of a lognormal and use the exact marginal variance log(1 + CV^2), so an
# interval here and a probability there rest on one distributional assumption.
#
# Plot with vanilla ggplot2; this file deliberately does not provide plot
# helpers.

# ===========================================================================
# generics
# ===========================================================================

#' Relative spawning stock trajectory
#'
#' @param object A run folder containing \code{Report.sso}, an \code{.rds} or
#'   \code{.Rdata} file holding an \code{\link[r4ss]{SS_output}} object, or such
#'   an object directly.
#' @param ... Additional arguments; see the methods.
#' @rdname ssfTrajectory
#' @export
setGeneric("ssfTrajectory",
  function(object, ...) standardGeneric("ssfTrajectory"))

#' Relative fishing mortality trajectory
#'
#' @param object A run folder containing \code{Report.sso}, an \code{.rds} or
#'   \code{.Rdata} file holding an \code{\link[r4ss]{SS_output}} object, or such
#'   an object directly.
#' @param ... Additional arguments; see the methods.
#' @rdname fTrajectory
#' @export
setGeneric("fTrajectory",
  function(object, ...) standardGeneric("fTrajectory"))

#' Lognormal confidence intervals
#'
#' @param object A trajectory \code{data.frame} or an \code{FLQuant} of point
#'   estimates.
#' @param ... Additional arguments; see the methods.
#' @rdname lognormalCI
#' @export
setGeneric("lognormalCI",
  function(object, ...) standardGeneric("lognormalCI"))

# ===========================================================================
# internal helpers
# ===========================================================================

#' Coerce assorted inputs to an SS_output object
#'
#' @param x Run folder, file path, or \code{\link[r4ss]{SS_output}} list.
#' @return An \code{SS_output}-shaped list.
#' @keywords internal
#' @noRd
ssLoad <- function(x) {

  isSS <- function(o)
    is.list(o) && all(c("timeseries", "derived_quants") %in% names(o))

  if (isSS(x)) return(x)

  if (is.character(x) && length(x) == 1 && dir.exists(x)) {
    if (!file.exists(file.path(x, "Report.sso")))
      stop("no Report.sso in ", x)
    hasCovar <- file.exists(file.path(x, "ss.cor")) ||
                file.exists(file.path(x, "ss3.cor"))
    if (!hasCovar)
      warning("no .cor file in ", x, ": standard errors will be unavailable, ",
              "so no intervals can be computed")
    return(r4ss::SS_output(x, verbose = FALSE, printstats = FALSE,
                           covar = hasCovar))
  }

  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    ext <- tolower(tools::file_ext(x))
    if (ext == "rds") {
      o <- readRDS(x)
      if (!isSS(o)) stop("the object in ", x, " is not SS_output-shaped")
      return(o)
    }
    if (ext %in% c("rdata", "rda")) {
      e   <- new.env()
      nms <- load(x, envir = e)
      for (nm in nms) {
        o <- get(nm, envir = e)
        if (isSS(o)) return(o)
      }
      stop("no SS_output-shaped object in ", x,
           " (loaded: ", paste(nms, collapse = ", "), ")")
    }
  }

  stop("unsupported input; supply a run folder, an .rds/.Rdata file, or an ",
       "SS_output object")
}

#' Pull labelled values and standard errors from derived_quants
#'
#' @param rep An \code{SS_output} object.
#' @param prefix Label prefix, e.g. \code{"SSB"}, \code{"Bratio"}, \code{"F"}.
#' @return A \code{data.frame} of \code{year}, \code{value}, \code{sd}, or
#'   \code{NULL} if no such labels exist.
#' @keywords internal
#' @noRd
ssSeries <- function(rep, prefix) {

  dq <- rep$derived_quants
  names(dq) <- tolower(names(dq))
  if (!"stddev" %in% names(dq)) dq$stddev <- NA_real_

  pat  <- sprintf("^%s_[0-9]{4}$", prefix)
  labs <- grep(pat, dq$label, value = TRUE)
  if (length(labs) == 0) return(NULL)

  yrs <- as.integer(sub(".*_", "", labs))
  o   <- order(yrs)
  i   <- match(labs[o], dq$label)

  data.frame(year  = yrs[o],
             value = as.numeric(dq$value[i]),
             sd    = as.numeric(dq$stddev[i]))
}

#' Single reference point from derived_quants, first label that is finite
#'
#' @param rep An \code{SS_output} object.
#' @param ... Candidate labels, tried in order.
#' @return Numeric, or \code{NA_real_}.
#' @keywords internal
#' @noRd
ssRef <- function(rep, ...) {

  dq <- rep$derived_quants
  names(dq) <- tolower(names(dq))

  for (lab in c(...)) {
    v <- as.numeric(dq$value[match(lab, dq$label)])
    if (length(v) == 1 && is.finite(v)) return(v)
  }
  NA_real_
}

# ===========================================================================
# trajectory methods
# ===========================================================================

#' @details
#' Values come from \code{derived_quants} rather than \code{timeseries} wherever
#' a label exists, because that is where the standard errors live. Pairing a
#' \code{timeseries} point estimate with a \code{derived_quants} standard error
#' is unsafe: the two need not refer to the same quantity once a model has
#' multiple seasons or areas, and the resulting CV is then wrong. The
#' \code{timeseries} table is a fallback only, in which case spawning biomass is
#' summed over areas within the spawning season rather than over all seasons, and
#' the \code{INIT} and \code{VIRG} equilibrium rows are dropped.
#'
#' Both \code{ssb} and \code{bratio} are returned, together with the standard
#' error of each. They are not interchangeable inputs to an interval: see
#' \code{\link{lognormalCI}}.
#'
#' @section Checking the depletion basis:
#' \code{Bratio} means \eqn{SSB/SSB_{MSY}} only when \code{starter.ss} sets the
#' depletion basis accordingly. The method reads \code{Bratio_label} and
#' \code{Bratio_denominator} from the \code{SS_output} object and warns if the
#' denominator is not \eqn{SSB_{MSY}}.
#'
#' @param label Prefix of the spawning-stock labels in \code{derived_quants}.
#'   Defaults to \code{"SSB"}; use \code{"SSF"} for models parameterised in
#'   spawning stock fecundity.
#' @return A \code{data.frame} with \code{year}, \code{value}, \code{valueSD},
#'   \code{ratio}, \code{ratioSD}, \code{refpt}, \code{quantity},
#'   \code{relative} and \code{source}.
#' @examples
#' \dontrun{
#' traj <- ssfTrajectory("path/to/run")
#' traj <- ssfTrajectory("HW2e_HW2d_InitFprior_1100t.Rdata")
#' }
#' @seealso \code{\link{fTrajectory}}, \code{\link{lognormalCI}},
#'   \code{\link{kobeProbs}}, \code{\link{getRuns}}
#' @rdname ssfTrajectory
#' @export
setMethod("ssfTrajectory", signature(object = "ANY"),
  function(object, label = "SSB") {

  rep <- ssLoad(object)

  denom <- rep$Bratio_denominator
  blab  <- rep$Bratio_label
  if (!is.null(denom) && !grepl("MSY", paste(denom, blab), ignore.case = TRUE))
    warning("the depletion basis is '", paste(na.omit(c(blab, denom)),
            collapse = " / "), "', not SSB_MSY. Bratio is therefore not ",
            "SSB/SSB_MSY and comparing it with one is not meaningful; set the ",
            "depletion basis in starter.ss to 2 if a Kobe ratio is wanted.")

  ssb    <- ssSeries(rep, label)
  bratio <- ssSeries(rep, "Bratio")

  if (is.null(ssb)) {

    warning("no ", label, "_<year> labels in derived_quants; falling back to ",
            "the timeseries table, where no standard errors are available")

    ts     <- rep$timeseries
    ssbCol <- if ("SpawnBio" %in% names(ts)) "SpawnBio" else "SSB"
    keep   <- !is.na(ts$Yr)
    if ("Era" %in% names(ts)) keep <- keep & !ts$Era %in% c("VIRG", "INIT")
    ts <- ts[keep, , drop = FALSE]

    if ("Seas" %in% names(ts) && length(unique(ts$Seas)) > 1) {
      bySeas <- tapply(ts[[ssbCol]], ts$Seas, sum, na.rm = TRUE)
      spawn  <- as.numeric(names(bySeas)[which.max(bySeas)])
      warning("multiple seasons: using season ", spawn, " as the spawning ",
              "season rather than summing over seasons")
      ts <- ts[ts$Seas == spawn, , drop = FALSE]
    }

    agg <- aggregate(list(value = ts[[ssbCol]]), by = list(year = ts$Yr),
                     FUN = sum, na.rm = TRUE)
    ssb <- data.frame(year = agg$year, value = agg$value, sd = NA_real_)
    src <- "timeseries"
  } else src <- "derived_quants"

  out <- data.frame(year = ssb$year, value = ssb$value, valueSD = ssb$sd,
                    quantity = "ssb", source = src, stringsAsFactors = FALSE)

  out$refpt <- ssRef(rep, paste0(label, "_MSY"), "SSB_MSY", "SSF_MSY")

  if (!is.null(bratio)) {
    i <- match(out$year, bratio$year)
    out$ratio    <- bratio$value[i]
    out$ratioSD  <- bratio$sd[i]
    out$relative <- "bratio"
  } else {
    out$ratio    <- out$value / out$refpt
    out$ratioSD  <- NA_real_
    out$relative <- "derived"
    warning("no Bratio_<year> labels: the ratio was formed by dividing by a ",
            "point estimate of the reference point, so it carries no standard ",
            "error of its own")
  }

  out[, c("year", "value", "valueSD", "ratio", "ratioSD", "refpt",
          "quantity", "relative", "source")]
})

#' @details
#' The relative fishing mortality series is taken from the \code{F_<year>}
#' labels in \code{derived_quants} together with their delta-method standard
#' errors. Whether \code{F_<year>} is \eqn{F/F_{MSY}} is read from
#' \code{F_std_basis}.
#'
#' @return A \code{data.frame} with the same columns as
#'   \code{\link{ssfTrajectory}}, with \code{quantity} equal to \code{"f"}.
#' @examples
#' \dontrun{
#' fTrajectory("path/to/run")
#' }
#' @seealso \code{\link{ssfTrajectory}}, \code{\link{lognormalCI}},
#'   \code{\link{kobeProbs}}, \code{\link{getRuns}}
#' @rdname fTrajectory
#' @export
setMethod("fTrajectory", signature(object = "ANY"), function(object) {

  rep <- ssLoad(object)

  f <- ssSeries(rep, "F")
  if (is.null(f))
    stop("no F_<year> labels in derived_quants; the sdreport may not cover ",
         "these years, or F reporting is switched off in starter.ss")

  basis  <- rep$F_std_basis
  isRel  <- !is.null(basis) &&
            grepl("fmsy|F_msy|/\\(F", paste(basis), ignore.case = TRUE)

  if (is.null(basis))
    warning("F_std_basis is absent from the SS_output object, so whether ",
            "F_<year> is already F/FMSY cannot be established; treating it as ",
            "absolute F. Check F_std_scaling in starter.ss.")

  out <- data.frame(year = f$year, value = f$value, valueSD = f$sd,
                    quantity = "f", source = "derived_quants",
                    stringsAsFactors = FALSE)

  out$refpt <- ssRef(rep, "annF_MSY", "Fstd_MSY")

  if (isRel) {
    out$ratio    <- f$value
    out$ratioSD  <- f$sd
    out$relative <- "reported"
  } else {
    out$ratio    <- f$value / out$refpt
    out$ratioSD  <- NA_real_
    out$relative <- "derived"
    warning("F_std_basis is '", paste(basis), "', so F_<year> is absolute F. ",
            "The ratio was formed by dividing by a point estimate of FMSY, ",
            "which ignores the uncertainty in FMSY and its correlation with F; ",
            "the resulting interval is conditional and too narrow. Setting ",
            "F_std_scaling = 2 in starter.ss gives F/FMSY with its own ",
            "standard error.")
  }

  out[, c("year", "value", "valueSD", "ratio", "ratioSD", "refpt",
          "quantity", "relative", "source")]
})

# ===========================================================================
# lognormalCI methods
# ===========================================================================

#' @details
#' The interval is \eqn{\hat{x}\exp(\pm z\sigma)} with
#' \eqn{\sigma = \sqrt{\log(1 + CV^2)}}, which treats the Stock Synthesis point
#' estimate as the median. Long form is returned, one row per year and level.
#'
#' @param levels Numeric vector of coverage probabilities, e.g.
#'   \code{c(0.5, 0.8, 0.95)}.
#' @param basis Either \code{"ratio"} or \code{"absolute"}.
#' @return \code{object} in long form with added columns \code{level},
#'   \code{lo}, \code{hi}, \code{ratioLo}, \code{ratioHi} and \code{basis}.
#' @examples
#' traj <- data.frame(year = c(2023, 2024), value = c(9000, 9500),
#'                    valueSD = c(1080, 1140), ratio = c(0.92, 0.95),
#'                    ratioSD = c(0.16, 0.166), refpt = 10000,
#'                    quantity = "ssb", relative = "bratio",
#'                    source = "derived_quants")
#' lognormalCI(traj, levels = 0.95)[2, c("ratioLo", "ratioHi")]
#' @seealso \code{\link{ssfTrajectory}}, \code{\link{fTrajectory}},
#'   \code{\link{kobeProbs}}, \code{\link{pAbove}}
#' @importFrom stats qnorm
#' @rdname lognormalCI
#' @export
setMethod("lognormalCI", signature(object = "data.frame"),
  function(object, levels = 0.95, basis = c("ratio", "absolute")) {

  basis <- match.arg(basis)

  if (any(!is.finite(levels) | levels <= 0 | levels >= 1))
    stop("'levels' must lie strictly between 0 and 1")

  need <- c("ratio", "ratioSD", "value", "valueSD", "refpt")
  if (!all(need %in% names(object)))
    stop("'object' must be trajectory output; missing: ",
         paste(setdiff(need, names(object)), collapse = ", "))

  if (basis == "ratio" && all(!is.finite(object$ratioSD))) {
    warning("no usable standard error on the ratio; falling back to ",
            "basis = 'absolute', which conditions on the reference point ",
            "being known and is generally too narrow")
    basis <- "absolute"
  }

  if (basis == "ratio") {
    x  <- object$ratio
    sd <- object$ratioSD
  } else {
    x  <- object$value
    sd <- object$valueSD
    if (any(is.finite(object$refpt)))
      warning("basis = 'absolute' ignores uncertainty in the reference point, ",
              "so the interval on the ratio is conditional and generally too ",
              "narrow; prefer basis = 'ratio'")
  }

  # Deterministic zero (e.g. F = 0 under zero catch): point mass at 0
  detZero <- is.finite(x) & x == 0 & is.finite(sd) & sd == 0
  bad <- !detZero & (!is.finite(sd) | sd <= 0 | !is.finite(x) | x <= 0)
  if (any(bad))
    warning(sum(bad), " years have a missing or zero standard error, or a ",
            "non-positive point estimate; their bounds are returned as NA ",
            "rather than as a zero-width interval")

  cv  <- ifelse(bad, NA_real_, ifelse(detZero, 0, sd / x))
  sig <- sqrt(log(1 + cv^2))

  build <- function(level) {
    z   <- stats::qnorm(1 - (1 - level) / 2)
    out <- object
    out$level <- level
    out$lo    <- x * exp(-z * sig)
    out$hi    <- x * exp( z * sig)
    if (basis == "ratio") {
      out$ratioLo <- out$lo
      out$ratioHi <- out$hi
    } else {
      out$ratioLo <- out$lo / object$refpt
      out$ratioHi <- out$hi / object$refpt
    }
    out$basis <- basis
    out
  }

  do.call(rbind, lapply(sort(levels), build))
})

#' @param sd An \code{FLQuant} of standard errors on the natural scale, with
#'   dimensions matching \code{object}.
#' @rdname lognormalCI
#' @export
setMethod("lognormalCI", signature(object = "FLQuant"),
  function(object, sd, levels = 0.95) {

  if (missing(sd))
    stop("the FLQuant method needs 'sd', an FLQuant of standard errors")
  if (!all(dim(object) == dim(sd)))
    stop("'object' and 'sd' must have the same dimensions")

  cv  <- sd / object
  sig <- sqrt(log(1 + cv^2))

  res <- list()
  for (level in sort(levels)) {
    z   <- stats::qnorm(1 - (1 - level) / 2)
    tag <- if (length(levels) == 1) c("lo", "hi")
           else paste0(c("lo", "hi"), round(level * 100))
    res[[tag[1]]] <- object * exp(-z * sig)
    res[[tag[2]]] <- object * exp( z * sig)
  }

  FLCore::FLQuants(res)
})
