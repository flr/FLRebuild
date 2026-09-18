# mvln.R - multivariate lognormal (MVLN) projections for Stock Synthesis
# FLRebuild
#
# Consolidated home for MVLN projection code:
#   * Monte Carlo draws (ssmvln and helpers) from the Walter–Winker / FLCandy
#     approximation, retained for path functionals and OM propagation
#   * Analytic orthant probabilities and Kobe II strategy matrices
#     (kobeGreen, kobeProbs, kobeMomentGrid, ...) — preferred for advice tables
#   * make* constructors for runs → moment grids → K2SM / rebuild references
#
# Prefer the analytic path for Kobe quadrant probabilities and K2SM tables;
# simulation only approximates a CDF that mvtnorm::pmvnorm evaluates directly.
#
# Trajectory extraction (ssfTrajectory, fTrajectory, lognormalCI) lives in
# trajectory.R; plot with vanilla ggplot2 rather than package plot helpers.
#
# ssmvln.R - multivariate lognormal approximation of SS3 Kobe posteriors
#
# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2
#
# ---------------------------------------------------------------------------
# MODIFIED WORK NOTICE (EUPL-1.2 Article 5, attribution right)
#
# This file is a modified version of ss-MVLN.R, moved from FLCandy to
# FLRebuild on 2026-09-17. Modifications on that date:
#   - roxygen2 documentation and NAMESPACE exports added
#   - cor2cov() supplied, having previously been an undeclared dependency
#   - shape and finiteness guards added, which stop or warn but do not alter
#     any returned value
#   - kobe_probs_mvln() retained as a deprecated alias of kobeProbsMVLN()
#
# No numerical behaviour was changed. The defects described under "Known
# issues" below are documented, not corrected, so that results computed with
# this file continue to match results computed with the FLCandy original.
#
# FLRebuild is distributed under GPL (>= 2), which is a Compatible Licence in
# the Appendix to the EUPL-1.2. Article 5 of the EUPL permits distribution of a
# derivative work combining EUPL and GPL material under the Compatible Licence,
# provided all copyright notices are kept intact and the modification is
# declared, as above.
# ---------------------------------------------------------------------------

#' Convert a correlation matrix to a covariance matrix
#'
#' Supplied here because the original \code{ss-MVLN.R} called \code{cor2cov()}
#' without declaring where it came from, so the file could not be sourced
#' standalone.
#'
#' @param cor A correlation matrix.
#' @param var A vector of variances, one per row of \code{cor}. Note the
#'   argument is variances, not standard deviations; the callers in this file
#'   pass \code{hat$cv}, which despite its name holds \eqn{\log(1 + CV^2)}.
#' @return A covariance matrix with the dimnames of \code{cor}.
#' @keywords internal
#' @noRd
cor2cov <- function(cor, var) {

  if (length(var) != nrow(cor))
    stop("length(var) must equal nrow(cor)")

  sd  <- sqrt(var)
  out <- outer(sd, sd) * cor
  dimnames(out) <- dimnames(cor)
  out
}

#' Extract point estimates and log-scale variances from derived_quants
#'
#' Selects the \code{Recr_}, \code{Bratio_} and \code{F_} rows of a Stock
#' Synthesis \code{derived_quants} table and returns their point estimates with
#' both natural- and log-scale dispersion.
#'
#' @section Column naming:
#' The returned column \code{cv} does \emph{not} hold a coefficient of
#' variation. It holds \eqn{\log(1 + CV^2)}, the log-scale variance, and
#' \code{stdLog} holds the squared CV. The names are inherited from the original
#' and are kept so that existing scripts continue to work. \code{cor2cov} is
#' called on \code{cv}, which is correct because \code{cor2cov} takes variances.
#'
#' The renaming step assumes \code{derived_quants} has exactly the columns that
#' \code{r4ss} produced in 2023, since it drops columns 4 and 5 by position and
#' then assigns five names. A shape check now stops with an informative error
#' rather than silently mislabelling columns if a newer \code{r4ss} returns a
#' different layout.
#'
#' @param hat A \code{derived_quants} table.
#' @return A \code{data.frame} with \code{label}, \code{hat}, \code{stdLog},
#'   \code{std} and \code{cv}.
#' @seealso \code{\link{ssmvln}}, \code{\link{kobeProbs}}
#' @export
sshat <- function(hat) {

  names(hat) <- tolower(names(hat))

  y <- rbind(hat[grep(paste("Recr",   "", sep = "_"), hat$label), ],
             hat[grep(paste("Bratio", "", sep = "_"), hat$label), ],
             hat[grep(paste("F",      "", sep = "_"), hat$label), ])

  y <- y[substr(y$label, 1, 2) %in% c("F_", "Br", "Re"), ]
  y <- subset(y, !label %in% c("Recr_Initial", "Recr_unfished",
                               "Recr_Virgin", "Ret_Catch_MSY"))

  y <- transform(y, cv2 = (y$stddev / y$value)^2,
                    var = log(1 + (y$stddev / y$value)^2))[, -c(4, 5)]

  if (ncol(y) != 5)
    stop("sshat expected 5 columns after dropping columns 4 and 5, got ",
         ncol(y), ". This selection is positional and depends on the r4ss ",
         "version; check the layout of derived_quants before proceeding.")

  names(y) <- c("label", "hat", "stdLog", "std", "cv")
  y
}

#' Correlation matrix of SS3 derived quantities
#'
#' Builds the correlation matrix among all \code{Recr_}, \code{Bratio_} and
#' \code{F_} quantities from a Stock Synthesis \code{CoVar} table.
#'
#' @section Known issue - absent pairs are set to zero:
#' Pairs with no \code{CoVar} entry are set to zero correlation. Zero is a
#' statement of independence, not of ignorance, so this silently substitutes a
#' different model for the one that was fitted. Two consequences matter. The
#' matrix may cease to be positive semi-definite, which
#' \code{\link{ssmvln}} then conceals by drawing with \code{method = "svd"}.
#' And a quantity wrongly given zero variance downstream produces a probability
#' of exactly 0 or 1, which is indistinguishable in a strategy matrix from a
#' genuine result. This behaviour is retained for compatibility but now warns.
#' \code{\link{kobeProbs}} takes the alternative approach of returning
#' \code{NA}.
#'
#' @param covar A \code{CoVar} table, coercible by
#'   \code{\link[data.table]{data.table}}.
#' @return A square correlation matrix with labels as dimnames.
#' @seealso \code{\link{ssmvln}}, \code{\link{kobeProbs}}
#' @importFrom data.table setnames dcast
#' @export
sscor <- function(covar) {

  data.table::setnames(covar, tolower(names(covar)))

  flag <- unique(sort(c(
    grep(paste("Recr",   "", sep = "_"), covar$label.i),
    grep(paste("Recr",   "", sep = "_"), covar$label.j),
    grep(paste("Bratio", "", sep = "_"), covar$label.i),
    grep(paste("Bratio", "", sep = "_"), covar$label.j),
    grep(paste("F",      "", sep = "_"), covar$label.i),
    grep(paste("F",      "", sep = "_"), covar$label.j))))

  cor <- covar[flag, c("label.i", "label.j", "corr")]

  flag <- substr(cor$label.i, 1, 2) %in% c("Re", "F_", "Br") &
          substr(cor$label.j, 1, 2) %in% c("Re", "F_", "Br")
  cor <- cor[flag, ]

  cor <- rbind(cor, transform(cor, label.i = label.j, label.j = label.i))
  cor <- subset(cor, !(is.na(label.i) | is.na(label.j)))

  cor  <- data.table::dcast(cor, label.i ~ label.j, value.var = "corr")
  dmns <- unname(unlist(cor[, 1]))

  cor <- as.matrix(cor[, -1])

  nMissing <- sum(is.na(cor)) - sum(is.na(diag(cor)))
  if (nMissing > 0)
    warning(nMissing, " off-diagonal correlations are absent from CoVar and ",
            "are being set to zero. Zero states independence, not ignorance, ",
            "and can render the matrix indefinite.")

  cor[is.na(cor)] <- 0
  diag(cor) <- 1

  dimnames(cor) <- list(dmns, dmns)
  cor
}

#' Multivariate lognormal draws from an SS3 covariance matrix
#'
#' Draws correlated time series of recruitment, relative biomass and relative
#' fishing mortality from the delta-method covariance matrix, as described by
#' Walter and Winker (2020). Serial correlation among years is preserved because
#' the whole matrix is used, not year-by-year marginals.
#'
#' @section When simulation is not needed:
#' For quadrant probabilities and Kobe II strategy matrices these draws are
#' unnecessary. The probability of a Kobe quadrant in a given year is a
#' bivariate normal orthant probability on the log scale and has a closed form;
#' see \code{\link{kobeProbs}} and \code{\link{kobeGreen}}. Simulating it only
#' introduces Monte Carlo error into a quantity that
#' \code{\link[mvtnorm]{pmvnorm}} evaluates directly. Use \code{ssmvln} when
#' joint draws across many years are genuinely wanted, for instance to propagate
#' trajectories into an operating model or to compute a functional of the whole
#' path.
#'
#' @section Known issue - the meaning of `new`:
#' \code{new = TRUE}, the default, passes \code{hat$cv} to \code{cor2cov}, which
#' is the log-scale variance \eqn{\log(1 + CV^2)} and is consistent with the
#' log-scale mean \code{log(hat$hat)}. \code{new = FALSE} passes
#' \code{hat$std^2}, a natural-scale variance, while the mean stays on the log
#' scale. The two are therefore not alternative parameterisations of one model:
#' \code{new = FALSE} is dimensionally inconsistent. Note that the original
#' \code{kobe_probs_mvln} defaulted to \code{new = FALSE} while \code{ssmvln}
#' defaulted to \code{TRUE}, so results from the two were not comparable. Both
#' defaults are preserved here; do not compare across them.
#'
#' @section Known issue - the covariance transform:
#' \code{cor2cov(cor, hat$cv)} sets the log-scale correlation equal to the
#' natural-scale correlation. The exact bivariate lognormal result is
#' \eqn{\rho_{\log} = \log(1 + \rho\, CV_u CV_v) / (\sigma_u \sigma_v)}. The two
#' agree closely at modest CVs and diverge as CVs grow. At large CV with strong
#' correlation the exact expression can imply \eqn{|\rho_{\log}| > 1}, meaning
#' the reported natural-scale correlation is unattainable for a bivariate
#' lognormal at all. \code{\link{kobeMomentGrid}} offers both transformations
#' and reports feasibility.
#'
#' @section Known issue - non-positive-definite matrices:
#' Drawing with \code{method = "svd"} succeeds on a matrix that is not positive
#' semi-definite, where a Cholesky factorisation would fail. Combined with the
#' zero-filling in \code{\link{sscor}} this means an inconsistent covariance
#' matrix produces draws without complaint. The default is retained for
#' reproducibility; a warning now reports the smallest eigenvalue when it is
#' negative.
#'
#' @param covar A \code{CoVar} table from \code{\link[r4ss]{SS_output}}.
#' @param hat A \code{derived_quants} table from the same object.
#' @param mc Number of draws. With \code{mc <= 1} the correlation matrix is
#'   returned instead of draws.
#' @param new Logical; see the note on \code{new} above. Defaults to
#'   \code{TRUE}, the log-scale variance.
#' @param method Passed to \code{\link[mvtnorm]{rmvnorm}}. Defaults to
#'   \code{"svd"} for compatibility.
#' @return A \code{data.table} in long form with columns \code{iter},
#'   \code{qname}, \code{year} and \code{data}; or a correlation matrix when
#'   \code{mc <= 1}.
#' @examples
#' \dontrun{
#' library(ss3diags)
#' data(sma)
#'
#' sscor(sma$CoVar)
#' sshat(sma$derived_quants)
#' ssmvln(sma$CoVar, sma$derived_quants)
#'
#' # the analytic equivalent for a single year, with no draws
#' kobeProbs(sma$derived_quants, sma$CoVar, 2017)
#' }
#' @references
#' Walter, J. and Winker, H. (2020) Projections to create Kobe 2 Strategy Matrix
#' using the multivariate log-normal approximation for Atlantic yellowfin tuna.
#' \emph{Collect. Vol. Sci. Pap. ICCAT} 76(6): 725-739.
#'
#' Kell, L., Rice, J., Courtney, D. and Winker, H. (2023) Approximation of Kobe
#' posteriors from Stock Synthesis for North Atlantic blue shark.
#' \emph{Collect. Vol. Sci. Pap. ICCAT} 80(4): 837-858.
#' @seealso \code{\link{kobeProbs}} for the analytic equivalent,
#'   \code{\link{kobeMomentGrid}} for the exact covariance transform
#' @importFrom mvtnorm rmvnorm
#' @importFrom data.table data.table melt tstrsplit setnames
#' @export
ssmvln <- function(covar, hat, mc = 5000, new = !FALSE, method = "svd") {

  names(hat)   <- tolower(names(hat))
  names(covar) <- tolower(names(covar))

  cor <- sscor(data.table::data.table(covar))

  if (mc <= 1) return(cor)

  hat <- data.table::data.table(hat)
  hat <- subset(hat, label %in% dimnames(cor)[[1]])
  hat <- sshat(hat)
  cor <- cor[hat$label, hat$label]

  if (new)
    cvr <- cor2cov(cor, hat$cv)
  else
    cvr <- cor2cov(cor, hat$std^2)

  nBad <- sum(!is.finite(cvr))
  if (nBad > 0)
    warning(nBad, " covariance entries are not finite and are being set to ",
            "zero. A zeroed diagonal entry gives that quantity zero variance, ",
            "and hence a probability of exactly 0 or 1 downstream.")

  cvr[!is.finite(cvr)] <- 0

  ev <- suppressWarnings(min(eigen(cvr, symmetric = TRUE,
                                   only.values = TRUE)$values))
  if (is.finite(ev) && ev < 0)
    warning("the covariance matrix is not positive semi-definite (smallest ",
            "eigenvalue ", signif(ev, 3), "); method = \"", method, "\" will ",
            "draw from it regardless, whereas a Cholesky factorisation would ",
            "fail")

  mvlmu <- log(hat$hat)
  names(mvlmu) <- dimnames(cor)[[1]]

  rtn <- exp(mvtnorm::rmvnorm(mc, mean = mvlmu, sigma = cvr, method = method))

  rtn <- data.table::data.table(rtn)
  names(rtn) <- dimnames(cor)[[1]]

  dat <- cbind(iter = seq(mc), rtn)
  dat <- data.table::melt(dat, id = "iter")
  dat <- dat[, c("variable", "year") :=
               data.table::tstrsplit(get("variable"), "_")]

  data.table::setnames(dat, c("variable", "value"), c("qname", "data"))
  dat
}

#' Covariance between stock and harvest as an FLQuant
#'
#' Extracts the diagonal of matched \code{Bratio_<year>} and \code{F_<year>}
#' pairs from a covariance matrix and returns it as an annual
#' \code{\link[FLCore]{FLQuant}}.
#'
#' @param covar A named covariance matrix, as returned by \code{\link{sscor}}.
#' @return An \code{FLQuant} of covariances by year.
#' @seealso \code{\link{ssmvln}}
#' @importFrom FLCore as.FLQuant
#' @export
covFLQ <- function(covar) {

  dt  <- cbind(expand.grid(dimnames(covar)), value = c(covar))
  dt  <- subset(dt, substr(ac(Var1), 8, nchar(ac(Var1))) ==
                    substr(ac(Var2), 3, nchar(ac(Var2))))
  cov <- cbind(dt, year = an(substr(ac(dt$Var2), 3, nchar(ac(dt$Var2)))))

  FLCore::as.FLQuant(data.frame(data = cov$value, year = cov$year))
}

#' Point estimates as FLQuants
#'
#' Splits a \code{derived_quants} table into recruitment, relative biomass and
#' relative fishing mortality series and returns them as
#' \code{\link[FLCore]{FLQuants}} named \code{recruits}, \code{stock} and
#' \code{harvest}.
#'
#' @param hat A \code{derived_quants} table.
#' @param data Name of the column to place in the \code{FLQuant}, e.g.
#'   \code{"hat"} for the point estimate or \code{"std"} for its standard error.
#' @return An \code{FLQuants} object.
#' @seealso \code{\link{sshat}}
#' @importFrom FLCore FLQuants
#' @export
hatFLQ <- function(hat, data = "hat") {

  r <- subset(hat, substr(label, 1, 5) == "Recr_")
  r <- transform(r, year = as.numeric(substr(label, 6, nchar(label))))

  b <- subset(hat, substr(label, 1, 7) == "Bratio_")
  b <- transform(b, year = as.numeric(substr(label, 8, nchar(label))))

  f <- subset(hat, substr(label, 1, 2) == "F_")
  f <- transform(f, year = as.numeric(substr(label, 3, nchar(label))))

  rtn <- rbind(cbind(qname = "recruits", r),
               cbind(qname = "stock",    b),
               cbind(qname = "harvest",  f))[, -2]
  rtn <- cbind(rtn[, c("qname", "year")], data = rtn[, data])

  as(rtn, "FLQuants")
}

#' MVLN draws as FLQuants
#'
#' Coerces the long-form output of \code{\link{ssmvln}} to
#' \code{\link[FLCore]{FLQuants}} named \code{recruits}, \code{stock} and
#' \code{harvest}.
#'
#' The names are assigned by position, so this assumes the coercion returns
#' exactly three quantities in that order. It will mislabel the result if
#' \code{\link{ssmvln}} was given a table containing other \code{F_} labels.
#'
#' @param x Long-form output of \code{\link{ssmvln}}.
#' @return An \code{FLQuants} object.
#' @seealso \code{\link{ssmvln}}
#' @export
mvlnFLQ <- function(x) {

  res <- as(x, "FLQuants")

  if (length(res) != 3)
    stop("expected 3 quantities to name recruits, stock and harvest, got ",
         length(res), "; the names would be assigned to the wrong series")

  names(res) <- c("recruits", "stock", "harvest")
  res
}

#' Kobe probabilities under the MVLN approximation
#'
#' Marginal and joint Kobe probabilities for one year, from the delta-method
#' covariance matrix.
#'
#' Despite its name this function draws no samples: it already evaluates the
#' joint probability with \code{\link[mvtnorm]{pmvnorm}}. It is therefore
#' superseded by \code{\link{kobeProbs}}, which does the same thing while also
#' returning all four quadrants, offering the exact lognormal covariance
#' transform, and returning \code{NA} rather than a degenerate 0 or 1 when a
#' standard error is missing.
#'
#' Note that \code{new} defaults to \code{FALSE} here but to \code{TRUE} in
#' \code{\link{ssmvln}}; see the note on \code{new} in that help page. The
#' inconsistency is preserved for backward compatibility.
#'
#' @param ss_out An \code{\link[r4ss]{SS_output}} object.
#' @param year Year for which probabilities are wanted. Defaults to the last
#'   year with a \code{Bratio_} label.
#' @param new Logical; see \code{\link{ssmvln}}. Defaults to \code{FALSE}.
#' @return A list with \code{year}, \code{p_F_lt_Fmsy}, \code{p_B_gt_Bmsy} and
#'   \code{p_good_quadrant}.
#' @examples
#' \dontrun{
#' kobeProbsMVLN(ss_out)
#' kobeProbsMVLN(ss_out, year = 2023)
#'
#' # preferred: all four quadrants, and the exact covariance transform
#' kobeProbs(ss_out$derived_quants, ss_out$CoVar, 2023)
#' }
#' @seealso \code{\link{kobeProbs}}, \code{\link{kobeGreen}}
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats pnorm
#' @export
kobeProbsMVLN <- function(ss_out, year = NULL, new = FALSE) {

  dq    <- ss_out$derived_quants
  covar <- ss_out$CoVar

  names(dq)    <- tolower(names(dq))
  names(covar) <- tolower(names(covar))

  if (is.null(year)) {
    brLabs <- grep("^Bratio_", dq$label, value = TRUE)
    year   <- max(as.numeric(sub("Bratio_", "", brLabs)))
  }

  corMat <- sscor(data.table::data.table(covar))
  hatDt  <- sshat(dq)

  hatDt  <- subset(hatDt, label %in% rownames(corMat))
  corMat <- corMat[hatDt$label, hatDt$label]

  if (new)
    cvr <- cor2cov(corMat, hatDt$cv)
  else
    cvr <- cor2cov(corMat, hatDt$std^2)

  cvr[!is.finite(cvr)] <- 0

  mvlmu <- log(hatDt$hat)
  names(mvlmu) <- hatDt$label

  labB <- sprintf("Bratio_%d", year)
  labF <- sprintf("F_%d",      year)

  if (!all(c(labB, labF) %in% names(mvlmu)))
    stop("requested year not found in Bratio_ or F_ labels")

  mu    <- mvlmu[c(labB, labF)]
  Sigma <- cvr[c(labB, labF), c(labB, labF), drop = FALSE]

  muB <- mu[1]; muF <- mu[2]
  sB2 <- Sigma[1, 1]; sF2 <- Sigma[2, 2]

  if (sB2 <= 0 || sF2 <= 0)
    warning("a zero variance for ", year, " will return a probability of ",
            "exactly 0 or 1; kobeProbs() returns NA in this situation")

  list(year            = year,
       p_F_lt_Fmsy     = stats::pnorm(0, mean = muF, sd = sqrt(sF2)),
       p_B_gt_Bmsy     = 1 - stats::pnorm(0, mean = muB, sd = sqrt(sB2)),
       p_good_quadrant = as.numeric(mvtnorm::pmvnorm(
         lower = c(0, -Inf), upper = c(Inf, 0),
         mean  = c(muB, muF), sigma = Sigma)))
}

#' @rdname kobeProbsMVLN
#' @details \code{kobe_probs_mvln} is a deprecated alias retained so that
#'   existing analysis scripts continue to run. Use \code{kobeProbsMVLN}, or
#'   preferably \code{\link{kobeProbs}}.
#' @export
kobe_probs_mvln <- function(ss_out, year = NULL, new = FALSE) {
  .Deprecated("kobeProbsMVLN")
  kobeProbsMVLN(ss_out, year = year, new = new)
}


# kobeAnalytic.R - analytic Kobe probabilities and rebuilding time
# FLRebuild
#
# Simulation-free evaluation of Kobe quadrant probabilities, the Kobe II
# strategy matrix, and rebuilding time, under the delta-MVLN approximation.

# ===========================================================================
# internal helpers: extraction from Stock Synthesis output
# ===========================================================================

#' Extract an MLE and standard error from derived_quants
#'
#' Missing, zero or non-finite standard errors return \code{NA} with a warning
#' rather than zero. A zero variance forces a degenerate probability of exactly
#' 0 or 1, which is indistinguishable from a genuine result once tabulated.
#'
#' @param dq \code{derived_quants} element of an \code{\link[r4ss]{SS_output}}
#'   object.
#' @param label Character, a \code{derived_quants} label, e.g. \code{"F_2024"}.
#' @return Named numeric vector with elements \code{value} and \code{se}.
#' @keywords internal
#' @noRd
kobeHat <- function(dq, label) {

  names(dq) <- tolower(names(dq))
  i <- match(label, dq$label)

  if (is.na(i)) {
    warning("label not found in derived_quants: ", label)
    return(c(value = NA_real_, se = NA_real_))
  }

  v  <- as.numeric(dq$value[i])
  se <- as.numeric(dq$stddev[i])

  if (!is.finite(v) || v <= 0) {
    warning("non-positive or non-finite value for ", label)
    return(c(value = NA_real_, se = NA_real_))
  }

  if (!is.finite(se) || se <= 0) {
    warning("missing or zero standard error for ", label,
            ": the log-scale variance would be zero, giving a probability of ",
            "exactly 0 or 1")
    return(c(value = v, se = NA_real_))
  }

  c(value = v, se = se)
}

#' Correlation between two derived quantities
#'
#' An absent pair returns \code{NA} with a warning. Setting an unknown
#' correlation to zero silently changes the model and can render the covariance
#' matrix indefinite.
#'
#' @param covar \code{CoVar} element of an \code{\link[r4ss]{SS_output}} object.
#' @param labelI,labelJ Character labels of the two quantities.
#' @return Numeric correlation, or \code{NA_real_}.
#' @keywords internal
#' @noRd
kobeRho <- function(covar, labelI, labelJ) {

  names(covar) <- tolower(names(covar))

  hit <- which((covar$label.i == labelI & covar$label.j == labelJ) |
               (covar$label.i == labelJ & covar$label.j == labelI))

  if (length(hit) == 0) {
    warning("no CoVar entry for ", labelI, " x ", labelJ,
            ": the correlation is unknown, not zero")
    return(NA_real_)
  }

  as.numeric(covar$corr[hit[1]])
}

#' Log-scale moments of relative biomass and fishing mortality
#'
#' Builds the bivariate normal mean vector and covariance matrix of
#' \eqn{\log(SSB/SSB_{MSY})} and \eqn{\log(F/F_{MSY})} for one year, from Stock
#' Synthesis delta-method output.
#'
#' Two transformations to the log scale are available. \code{"exact"} uses the
#' lognormal result \eqn{cov = \log(1 + \rho\,CV_u CV_v)}. \code{"current"}
#' reproduces the behaviour of \code{ssmvln()}, which sets the log-scale
#' correlation equal to the natural-scale correlation. The two agree closely at
#' modest CVs and diverge as CVs grow; at large CV combined with strong
#' correlation the exact expression can imply \eqn{|\rho_{\log}| > 1}, meaning
#' the reported natural-scale correlation is not attainable for a bivariate
#' lognormal. That case is reported through \code{feasible}.
#'
#' The SS3 point estimate is treated as the median of the lognormal, so the
#' log-scale mean is \code{log(value)}. The mean of the lognormal therefore
#' exceeds the MLE.
#'
#' @param dq \code{derived_quants} element of an \code{\link[r4ss]{SS_output}}
#'   object.
#' @param covar \code{CoVar} element of the same object.
#' @param year Numeric year.
#' @param transform Either \code{"exact"} or \code{"current"}.
#' @param labB,labF Character labels, by default \code{Bratio_<year>} and
#'   \code{F_<year>}. Requires \code{starter.ss} to set the depletion basis to
#'   \eqn{SSB_{MSY}} and \code{F_std_scaling} to \eqn{F/F_{MSY}}; otherwise the
#'   labels do not carry the assumed meaning.
#' @return List with \code{mu}, \code{Sigma}, \code{sigmaB}, \code{sigmaF},
#'   \code{rhoLog}, \code{feasible}, \code{year}, \code{bratio} and
#'   \code{fratio}.
#' @keywords internal
#' @noRd
kobeMoments <- function(dq, covar, year,
                        transform = c("exact", "current"),
                        labB = sprintf("Bratio_%d", year),
                        labF = sprintf("F_%d", year)) {

  transform <- match.arg(transform)

  b   <- kobeHat(dq, labB)
  f   <- kobeHat(dq, labF)
  rho <- kobeRho(covar, labB, labF)

  empty <- list(mu = c(NA_real_, NA_real_), Sigma = matrix(NA_real_, 2, 2),
                sigmaB = NA_real_, sigmaF = NA_real_, rhoLog = NA_real_,
                feasible = NA, year = year, bratio = NA_real_,
                fratio = NA_real_)

  if (any(is.na(c(b, f, rho)))) return(empty)

  cvB <- b["se"] / b["value"]
  cvF <- f["se"] / f["value"]

  s2B <- log(1 + cvB^2)
  s2F <- log(1 + cvF^2)
  sB  <- sqrt(s2B)
  sF  <- sqrt(s2F)

  covLog <- switch(transform,
    exact   = log(1 + rho * cvB * cvF),
    current = rho * sB * sF)

  rhoLog   <- covLog / (sB * sF)
  feasible <- is.finite(rhoLog) && abs(rhoLog) < 1

  if (!feasible)
    warning("year ", year, ": implied log-scale correlation ",
            round(rhoLog, 3), " lies outside (-1, 1); the reported ",
            "natural-scale correlation and CVs are not jointly attainable ",
            "for a bivariate lognormal")

  list(mu = c(unname(log(b["value"])), unname(log(f["value"]))),
       Sigma = matrix(c(s2B, covLog, covLog, s2F), 2, 2),
       sigmaB = sB, sigmaF = sF, rhoLog = rhoLog, feasible = feasible,
       year = year, bratio = unname(b["value"]), fratio = unname(f["value"]))
}

# ===========================================================================
# exported functions
# ===========================================================================

#' Kobe quadrant probability from log-scale moments
#'
#' Evaluates the probability that a stock lies in a Kobe quadrant as an orthant
#' probability of a bivariate normal on the log scale. Under the delta-MVLN
#' approximation of Walter and Winker (2020) this quantity has a closed form, so
#' no Monte Carlo simulation is required: simulation only approximates a
#' cumulative distribution function that \code{\link[mvtnorm]{pmvnorm}}
#' evaluates directly.
#'
#' For the green quadrant,
#' \deqn{P = \Phi_2\left(\frac{\mu_B}{\sigma_B},\,
#'   -\frac{\mu_F}{\sigma_F};\, -\rho\right)}
#' with \eqn{\Phi_2} the standard bivariate normal distribution function.
#'
#' @param muB,muF Means of \eqn{\log(SSB/SSB_{MSY})} and
#'   \eqn{\log(F/F_{MSY})}.
#' @param sigmaB,sigmaF Standard deviations on the log scale.
#' @param rho Correlation on the log scale.
#' @param quadrant One of \code{"green"} (not overfished, no overfishing),
#'   \code{"red"} (overfished and overfishing), \code{"yellow"} (overfished, no
#'   overfishing) or \code{"orange"} (not overfished, overfishing).
#' @return Numeric probability, or \code{NA_real_} if any moment is not finite
#'   or the correlation is not inside \eqn{(-1, 1)}.
#' @examples
#' # 2026 North Atlantic shortfin mako terminal year, moments implied by the
#' # reported point probabilities
#' kobeGreen(muB = log(0.95), muF = log(0.72),
#'           sigmaB = 0.174, sigmaF = 0.384, rho = -0.652)
#'
#' # the four quadrants sum to one
#' sum(sapply(c("green", "red", "yellow", "orange"), function(q)
#'   kobeGreen(log(0.95), log(0.72), 0.174, 0.384, -0.652, quadrant = q)))
#' @references
#' Walter, J. and Winker, H. (2020) Projections to create Kobe 2 Strategy Matrix
#' using the multivariate log-normal approximation for Atlantic yellowfin tuna.
#' \emph{Collect. Vol. Sci. Pap. ICCAT} 76(6): 725-739.
#'
#' Kell, L., Rice, J., Courtney, D. and Winker, H. (2023) Approximation of Kobe
#' posteriors from Stock Synthesis for North Atlantic blue shark.
#' \emph{Collect. Vol. Sci. Pap. ICCAT} 80(4): 837-858.
#' @seealso \code{\link{kobeProbs}}, \code{\link{kobeRebuildTime}}
#' @importFrom mvtnorm pmvnorm
#' @export
kobeGreen <- function(muB, muF, sigmaB, sigmaF, rho,
                      quadrant = c("green", "red", "yellow", "orange")) {

  quadrant <- match.arg(quadrant)

  if (any(!is.finite(c(muB, muF, sigmaB, sigmaF, rho)))) return(NA_real_)
  if (abs(rho) >= 1 || sigmaB <= 0 || sigmaF <= 0) return(NA_real_)

  Sigma <- matrix(c(sigmaB^2, rho * sigmaB * sigmaF,
                    rho * sigmaB * sigmaF, sigmaF^2), 2, 2)

  bounds <- switch(quadrant,
    green  = list(lower = c(0, -Inf),    upper = c(Inf, 0)),
    red    = list(lower = c(-Inf, 0),    upper = c(0, Inf)),
    yellow = list(lower = c(-Inf, -Inf), upper = c(0, 0)),
    orange = list(lower = c(0, 0),       upper = c(Inf, Inf)))

  as.numeric(mvtnorm::pmvnorm(lower = bounds$lower, upper = bounds$upper,
                              mean = c(muB, muF), sigma = Sigma))
}

#' Kobe quadrant probabilities for one year
#'
#' All four quadrant probabilities and both marginals for a single year,
#' evaluated analytically from Stock Synthesis delta-method output.
#'
#' The returned \code{sigmaLogF} is a diagnostic as well as a parameter. If it
#' falls markedly below its value in the terminal assessment year, the
#' projection is treating fishing mortality as better determined than the
#' assessment estimates it to be, which drives
#' \eqn{P(F \le F_{MSY})} towards one. See \code{\link{kobeMomentGrid}}.
#'
#' @param dq \code{derived_quants} element of an \code{\link[r4ss]{SS_output}}
#'   object.
#' @param covar \code{CoVar} element of the same object.
#' @param year Numeric year.
#' @param transform Either \code{"exact"} or \code{"current"}; see Details of
#'   \code{\link{kobeMomentGrid}}.
#' @return A one-row \code{data.frame} with \code{year}, \code{sigmaLogB},
#'   \code{sigmaLogF}, \code{rhoLog}, \code{pNotOverfished},
#'   \code{pNoOverfishing}, and \code{green}, \code{red}, \code{yellow} and
#'   \code{orange}.
#' @examples
#' \dontrun{
#' ss <- r4ss::SS_output("path/to/run")
#' kobeProbs(ss$derived_quants, ss$CoVar, 2024)
#' }
#' @seealso \code{\link{kobeGreen}}, \code{\link{kobeMomentGrid}}
#' @importFrom stats pnorm
#' @export
kobeProbs <- function(dq, covar, year, transform = c("exact", "current")) {

  transform <- match.arg(transform)
  m <- kobeMoments(dq, covar, year, transform)

  empty <- data.frame(year = year, sigmaLogB = NA_real_, sigmaLogF = NA_real_,
                      rhoLog = NA_real_, pNotOverfished = NA_real_,
                      pNoOverfishing = NA_real_, green = NA_real_,
                      red = NA_real_, yellow = NA_real_, orange = NA_real_)

  if (any(is.na(m$mu))) return(empty)

  q <- function(w) kobeGreen(m$mu[1], m$mu[2], m$sigmaB, m$sigmaF, m$rhoLog, w)

  data.frame(year = year,
             sigmaLogB = m$sigmaB, sigmaLogF = m$sigmaF, rhoLog = m$rhoLog,
             pNotOverfished = 1 - stats::pnorm(0, m$mu[1], m$sigmaB),
             pNoOverfishing = stats::pnorm(0, m$mu[2], m$sigmaF),
             green = q("green"), red = q("red"),
             yellow = q("yellow"), orange = q("orange"))
}

#' Log-scale moments across catch levels and years
#'
#' Assembles the bivariate normal moments of \eqn{\log(SSB/SSB_{MSY})} and
#' \eqn{\log(F/F_{MSY})} for every combination of projection run and year. This
#' is the object interpolated by \code{\link{kobeInterp}} and consumed by
#' \code{\link{kobeRebuildCurve}}.
#'
#' Because \code{starter.ss} can be set to extend the sdreport to
#' \code{endyr + Nforecastyrs}, standard errors are available for every forecast
#' year, so moments need not be interpolated in time. Supply the full annual
#' sequence rather than the years tabulated in a strategy matrix.
#'
#' All runs must share an identical \code{forecast.ss} apart from the catch
#' vector. Forecast settings can alter the between-year correlations, so
#' interpolating across runs that differ in other respects mixes distinct
#' models.
#'
#' Run x year combinations with missing moments are retained as \code{NA} and
#' reported through a warning. The usual cause is a zero or absent standard
#' error on a forecast-year \code{F_} label, which would otherwise produce a
#' probability of exactly one.
#'
#' @param runs Named list of \code{\link[r4ss]{SS_output}} objects, one per
#'   catch level. Names must be the catch levels and coercible to numeric.
#' @param years Numeric vector of years.
#' @param transform Either \code{"exact"}, using the lognormal result
#'   \eqn{\log(1 + \rho\,CV_u CV_v)}, or \code{"current"}, reproducing
#'   \code{ssmvln()} by setting the log-scale correlation equal to the
#'   natural-scale one.
#' @return A \code{data.frame} with \code{catch}, \code{year}, \code{muB},
#'   \code{muF}, \code{sigmaB}, \code{sigmaF}, \code{rhoLog} and
#'   \code{feasible}.
#' @examples
#' \dontrun{
#' runs <- list("0" = ss0, "250" = ss250, "1000" = ss1000, "1500" = ss1500)
#' grid <- kobeMomentGrid(runs, 2025:2070)
#'
#' # is projected F uncertainty consistent with the assessment?
#' with(subset(grid, catch == 250),
#'      plot(year, sigmaF / sigmaF[year == 2024], type = "l"))
#' }
#' @seealso \code{\link{kobeInterp}}, \code{\link{kobeRebuildCurve}}
#' @importFrom stats complete.cases
#' @export
kobeMomentGrid <- function(runs, years, transform = c("exact", "current")) {

  transform <- match.arg(transform)

  catch <- suppressWarnings(as.numeric(names(runs)))
  if (length(catch) == 0 || any(is.na(catch)))
    stop("names(runs) must be the catch levels, coercible to numeric")

  out <- do.call(rbind, lapply(seq_along(runs), function(k) {
    ss <- runs[[k]]
    do.call(rbind, lapply(years, function(y) {
      m <- kobeMoments(ss$derived_quants, ss$CoVar, y, transform)
      data.frame(catch = catch[k], year = y,
                 muB = m$mu[1], muF = m$mu[2],
                 sigmaB = m$sigmaB, sigmaF = m$sigmaF,
                 rhoLog = m$rhoLog, feasible = m$feasible)
    }))
  }))

  out <- out[order(out$catch, out$year), ]
  rownames(out) <- NULL

  bad <- !stats::complete.cases(
    out[, c("muB", "muF", "sigmaB", "sigmaF", "rhoLog")])

  if (any(bad))
    warning(sum(bad), " run x year combinations have missing moments. ",
            "Check for zero or absent standard errors on forecast-year ",
            "F_ labels.")

  out
}

#' Interpolate Kobe moments across catch levels
#'
#' Interpolates the log-scale moments across catch and then evaluates the
#' quadrant probability exactly, rather than interpolating the probabilities
#' themselves.
#'
#' The probability is a bounded, sigmoid-shaped function of the moments, so
#' interpolating it inherits that curvature as error. The mean
#' \eqn{\mu_B = \log(SSB/SSB_{MSY})} is close to linear in cumulative removals,
#' and is therefore the better quantity to interpolate. Evaluating
#' \code{\link{kobeGreen}} at interpolated moments also guarantees that the
#' result remains a probability from a valid distribution.
#'
#' Monotone Hermite interpolation (Fritsch and Carlson, 1980, via
#' \code{\link[stats]{splinefun}} method \code{"monoH.FC"}) is used in
#' preference to an unconstrained spline: the probability is monotone decreasing
#' in catch and the surface is nearly flat at low catch, where a natural cubic
#' will overshoot. With two runs the interpolation is linear.
#'
#' Extrapolation beyond the fitted range is warned about rather than refused.
#' The catch-probability surface can have a kink where the forecast hits an F
#' cap or the stock approaches collapse, and interpolation across such a kink is
#' not valid.
#'
#' @param grid Output of \code{\link{kobeMomentGrid}}.
#' @param catch Numeric vector of catch levels at which to evaluate.
#' @param years Numeric vector of years; defaults to all years in \code{grid}.
#' @param quadrant Kobe quadrant, passed to \code{\link{kobeGreen}}.
#' @return A \code{data.frame} with \code{catch}, \code{year}, the interpolated
#'   moments, and \code{green} holding the probability of the requested
#'   quadrant.
#' @examples
#' \dontrun{
#' grid <- kobeMomentGrid(runs, 2025:2070)
#' kobeInterp(grid, catch = c(1100, 1174), years = c(2040, 2050, 2070))
#' }
#' @references
#' Fritsch, F.N. and Carlson, R.E. (1980) Monotone piecewise cubic
#' interpolation. \emph{SIAM Journal on Numerical Analysis} 17(2): 238-246.
#' @seealso \code{\link{kobeMomentGrid}}, \code{\link{kobeRebuildCurve}}
#' @importFrom stats approx splinefun
#' @export
kobeInterp <- function(grid, catch, years = NULL,
                       quadrant = c("green", "red", "yellow", "orange")) {

  quadrant <- match.arg(quadrant)
  if (is.null(years)) years <- sort(unique(grid$year))

  interp1 <- function(x, y, xout) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    o <- order(x); x <- x[o]; y <- y[o]
    keep <- !duplicated(x); x <- x[keep]; y <- y[keep]
    if (length(x) < 2)  return(rep(NA_real_, length(xout)))
    if (length(x) == 2) return(stats::approx(x, y, xout = xout, rule = 2)$y)
    stats::splinefun(x, y, method = "monoH.FC")(xout)
  }

  rng <- range(grid$catch[is.finite(grid$muB)], na.rm = TRUE)
  if (any(catch < rng[1] | catch > rng[2]))
    warning("extrapolating beyond the fitted catch range (", rng[1], "-",
            rng[2], "); the catch-probability surface can have a kink where ",
            "the forecast hits an F cap or the stock approaches collapse")

  res <- do.call(rbind, lapply(years, function(y) {
    g  <- grid[grid$year == y, ]
    mb <- interp1(g$catch, g$muB,    catch)
    mf <- interp1(g$catch, g$muF,    catch)
    sb <- interp1(g$catch, g$sigmaB, catch)
    sf <- interp1(g$catch, g$sigmaF, catch)
    rl <- interp1(g$catch, g$rhoLog, catch)
    data.frame(catch = catch, year = y,
               muB = mb, muF = mf, sigmaB = sb, sigmaF = sf, rhoLog = rl,
               green = mapply(kobeGreen, mb, mf, sb, sf, rl,
                              MoreArgs = list(quadrant = quadrant)))
  }))

  res[order(res$catch, res$year), ]
}

#' Rebuilding time from an annual probability series
#'
#' Rebuilding time is defined as the first year from which the probability
#' \emph{remains} at or above the target, that is the last crossing of the
#' target rather than the first.
#'
#' A first-crossing rule is unsafe because the probability need not be monotone
#' in time. For an age-structured stock recovering from a truncated age
#' distribution the probability typically falls to a minimum some years into the
#' projection before rising, at every catch level including zero. A first
#' crossing can therefore be picked up before the transient rather than after
#' it, reporting as rebuilt a stock that is on its way down.
#'
#' The returned \code{dPdt} is a conditioning diagnostic. Rebuilding time is
#' ill-determined as \code{dPdt} approaches zero, which happens in two places:
#' inside the transient, where the derivative changes sign, and on the long flat
#' approach to the equilibrium, where a small change in probability moves the
#' crossing year by decades. In both regions the probability at fixed years is
#' the more defensible summary.
#'
#' @param prob A \code{data.frame} with columns \code{year} and \code{green},
#'   such as one catch level of \code{\link{kobeInterp}} output.
#' @param target Target probability, on \eqn{(0, 1)}. Recommendation 21-09 for
#'   North Atlantic shortfin mako specifies 0.60 to 0.70.
#' @return A one-row \code{data.frame} with \code{target}, \code{rebuildYear},
#'   \code{nCrossings}, \code{dPdt} in units of probability per year, and
#'   \code{attained}, which is \code{FALSE} when the target is not held at the
#'   end of the series.
#' @examples
#' # a probability series with a transient dip, as seen in age-structured
#' # rebuilding projections
#' prob <- data.frame(year  = c(2027, 2030, 2035, 2040, 2045, 2050, 2070),
#'                    green = c(0.50, 0.44, 0.55, 0.80, 0.92, 0.97, 1.00))
#' kobeRebuildTime(prob, target = 0.60)
#'
#' # at a lower target the series crosses twice; the last crossing is taken
#' kobeRebuildTime(prob, target = 0.50)
#' @seealso \code{\link{kobeRebuildCurve}}, \code{\link{kobeExchangeRate}}
#' @export
kobeRebuildTime <- function(prob, target = 0.6) {

  if (!all(c("year", "green") %in% names(prob)))
    stop("'prob' must have columns 'year' and 'green'")
  if (!is.finite(target) || target <= 0 || target >= 1)
    stop("'target' must lie strictly between 0 and 1")

  prob <- prob[order(prob$year), ]
  y <- prob$year
  p <- prob$green
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]; p <- p[ok]

  empty <- function(n, att)
    data.frame(target = target, rebuildYear = NA_real_, nCrossings = n,
               dPdt = NA_real_, attained = att)

  if (length(y) < 2) return(empty(NA_integer_, NA))

  cross <- numeric(0)
  for (i in seq_len(length(y) - 1)) {
    a <- p[i] - target
    b <- p[i + 1] - target
    if (a * b <= 0 && p[i] != p[i + 1])
      cross <- c(cross,
                 y[i] + (y[i + 1] - y[i]) * (target - p[i]) / (p[i + 1] - p[i]))
  }

  if (p[length(p)] < target) return(empty(length(cross), FALSE))

  # at or above the target throughout the series
  if (length(cross) == 0)
    return(data.frame(target = target, rebuildYear = min(y), nCrossings = 0L,
                      dPdt = NA_real_, attained = TRUE))

  tStar <- max(cross)
  i <- min(max(which(y <= tStar)), length(y) - 1)
  dPdt <- (p[i + 1] - p[i]) / (y[i + 1] - y[i])

  if (abs(dPdt) < 0.005)
    warning("rebuilding time is poorly determined: the probability changes by ",
            round(dPdt * 100, 2), " percentage points per year at the ",
            "crossing; report the probability at fixed years instead")

  data.frame(target = target, rebuildYear = tStar,
             nCrossings = length(cross), dPdt = dPdt, attained = TRUE)
}

#' Rebuilding time as a function of catch
#'
#' Applies \code{\link{kobeRebuildTime}} to each catch level of a moment grid,
#' interpolating the moments where a catch level was not run.
#'
#' @param grid Output of \code{\link{kobeMomentGrid}}.
#' @param catch Numeric vector of catch levels.
#' @param years Numeric vector of years; defaults to all years in \code{grid}.
#'   Supply the full annual sequence, not the years of a strategy matrix, since
#'   coarse year spacing determines the resolution of the crossing.
#' @param target Target probability.
#' @return A \code{data.frame} with \code{catch} and the columns returned by
#'   \code{\link{kobeRebuildTime}}.
#' @examples
#' \dontrun{
#' grid  <- kobeMomentGrid(runs, 2025:2070)
#' curve <- kobeRebuildCurve(grid, seq(0, 1800, by = 25), target = 0.60)
#' plot(curve$catch, curve$rebuildYear, type = "l",
#'      xlab = "Total removals (t)", ylab = "Year P(green) >= 0.60")
#' }
#' @seealso \code{\link{kobeRebuildTime}}, \code{\link{kobeExchangeRate}}
#' @export
kobeRebuildCurve <- function(grid, catch, years = NULL, target = 0.6) {

  surf <- kobeInterp(grid, catch, years)

  out <- do.call(rbind, lapply(catch, function(c0)
    cbind(catch = c0, kobeRebuildTime(surf[surf$catch == c0, ], target))))

  rownames(out) <- NULL
  out
}

#' Exchange rate between catch and rebuilding time
#'
#' The derivative of rebuilding time with respect to catch, which expresses how
#' many years of delay are bought by an increment of catch. This is the quantity
#' of interest for risk-equivalent advice, since it converts a change in a catch
#' limit into a change in the time taken to meet a probability criterion.
#'
#' Two estimates are returned and should agree. \code{dTdC} is a central
#' difference on the rebuilding curve. \code{dTdCimplicit} applies the implicit
#' function theorem at the crossing,
#' \deqn{\frac{dT}{dC} = -\frac{\partial P/\partial C}{\partial P/\partial t}}
#' which makes the source of any instability explicit: the derivative diverges
#' as \eqn{\partial P/\partial t} approaches zero. Disagreement between the two
#' indicates that the crossing is poorly conditioned, or that the interpolated
#' surface is not smooth over the step used.
#'
#' @param grid Output of \code{\link{kobeMomentGrid}}.
#' @param catch Numeric vector of catch levels at which to differentiate.
#' @param years Numeric vector of years; defaults to all years in \code{grid}.
#' @param target Target probability.
#' @param h Step in catch units for the central difference. Defaults to
#'   1/200 of the fitted catch range.
#' @return A \code{data.frame} with \code{catch}, \code{rebuildYear},
#'   \code{dPdt}, \code{dTdC}, \code{dTdCimplicit} and \code{yearsPer100t}.
#' @examples
#' \dontrun{
#' grid <- kobeMomentGrid(runs, 2025:2070)
#' kobeExchangeRate(grid, catch = c(250, 1000, 1250, 1500), target = 0.60)
#' }
#' @seealso \code{\link{kobeRebuildCurve}}, \code{\link{kobeTargetCatch}}
#' @export
kobeExchangeRate <- function(grid, catch, years = NULL, target = 0.6,
                             h = NULL) {

  if (is.null(h)) h <- max(1, diff(range(grid$catch, na.rm = TRUE)) / 200)

  out <- do.call(rbind, lapply(catch, function(c0) {

    t0 <- kobeRebuildTime(kobeInterp(grid, c0,     years), target)
    tm <- kobeRebuildTime(kobeInterp(grid, c0 - h, years), target)
    tp <- kobeRebuildTime(kobeInterp(grid, c0 + h, years), target)

    fd <- if (is.finite(tm$rebuildYear) && is.finite(tp$rebuildYear))
            (tp$rebuildYear - tm$rebuildYear) / (2 * h) else NA_real_

    imp <- NA_real_
    if (is.finite(t0$rebuildYear) && is.finite(t0$dPdt) && t0$dPdt != 0) {
      yr <- round(t0$rebuildYear)
      sm <- kobeInterp(grid, c0 - h, yr)$green
      sp <- kobeInterp(grid, c0 + h, yr)$green
      if (is.finite(sm) && is.finite(sp))
        imp <- -((sp - sm) / (2 * h)) / t0$dPdt
    }

    data.frame(catch = c0, rebuildYear = t0$rebuildYear, dPdt = t0$dPdt,
               dTdC = fd, dTdCimplicit = imp, yearsPer100t = fd * 100)
  }))

  rownames(out) <- NULL
  out
}

#' Catch achieving a target probability by a given year
#'
#' Inverts the interpolated probability surface to find the catch level at which
#' a stated probability is met by a stated year. This is the risk-equivalence
#' counterpart of \code{\link{kobeRebuildCurve}}: instead of asking when a catch
#' rebuilds the stock, it asks what catch rebuilds it by a chosen date.
#'
#' @param grid Output of \code{\link{kobeMomentGrid}}.
#' @param year Numeric year by which the target must be met.
#' @param target Target probability.
#' @param range Numeric vector of length two bounding the search in catch units;
#'   defaults to the fitted catch range.
#' @return Numeric catch level. \code{NA_real_} if the target is not met even at
#'   the lower bound of \code{range}; the upper bound if it is met throughout.
#' @examples
#' \dontrun{
#' grid <- kobeMomentGrid(runs, 2025:2070)
#' sapply(c(2040, 2050, 2070), function(y)
#'   kobeTargetCatch(grid, y, target = 0.60))
#' }
#' @seealso \code{\link{kobeExchangeRate}}
#' @importFrom stats uniroot
#' @export
kobeTargetCatch <- function(grid, year, target = 0.6, range = NULL) {

  if (is.null(range))
    range <- range(grid$catch[is.finite(grid$muB)], na.rm = TRUE)

  f <- function(c0) kobeInterp(grid, c0, year)$green - target

  lo <- f(range[1])
  hi <- f(range[2])

  if (!is.finite(lo) || !is.finite(hi)) return(NA_real_)
  if (lo < 0) return(NA_real_)
  if (hi > 0) return(range[2])

  stats::uniroot(f, range)$root
}

# ===========================================================================
# make* constructors (sma-rebuild workflow as S4 generics / methods)
# ===========================================================================

#' Marginal probability that a ratio exceeds a threshold
#'
#' Under the same lognormal assumption used by \code{\link{lognormalCI}} and
#' the MVLN projections,
#' \eqn{P(x > t) = \Phi(\log(x / t) / \sigma)} with
#' \eqn{\sigma = \sqrt{\log(1 + CV^2)}}.
#'
#' @param ratio Point estimate of the ratio (e.g. \eqn{SSF/SSF_{MSY}}).
#' @param ratioSD Standard error on the natural scale.
#' @param threshold Threshold on the ratio scale; defaults to 1.
#' @return Numeric probability, vectorised over \code{ratio} and \code{ratioSD}.
#' @examples
#' # 2026 North Atlantic shortfin mako terminal year
#' pAbove(0.95, 0.16656)   # P(SSF > SSFMSY) ≈ 0.384
#' pAbove(0.72, 0.2836)    # P(F > FMSY) ≈ 0.196
#' @seealso \code{\link{kobeProbs}}, \code{\link{lognormalCI}}
#' @importFrom stats pnorm
#' @export
pAbove <- function(ratio, ratioSD, threshold = 1) {
  cv    <- ratioSD / ratio
  sigma <- sqrt(log(1 + cv^2))
  stats::pnorm(log(ratio / threshold) / sigma)
}

#' Load Stock Synthesis projection runs
#'
#' @param object Named character vector of run folders or \code{.rds}/\code{.Rdata}
#'   paths, a named list of \code{SS_output} objects, or a directory whose
#'   \code{*.Rdata} files are named \code{*_\\{catch\\}t.Rdata}.
#' @param ... Unused.
#' @return A named list of \code{SS_output}-shaped objects; names are catch
#'   levels coercible to numeric (required by \code{\link{kobeMomentGrid}}).
#' @examples
#' \dontrun{
#' runs <- getRuns(c(
#'   "0" = "path/to/HW2e_..._0t.Rdata",
#'   "250" = "path/to/HW2e_..._250t.Rdata"))
#' runs <- getRuns("P:/.../SSoutput")  # discovers *_\\{catch\\}t.Rdata
#' }
#' @seealso \code{\link{kobeMomentGrid}}, \code{\link{makeK2SM}}
#' @export
setGeneric("getRuns", function(object, ...) standardGeneric("getRuns"))

#' @rdname getRuns
#' @export
setMethod("getRuns", signature(object = "character"),
  function(object, pattern = "_([0-9]+)t\\.Rdata$", ...) {

  if (length(object) == 1L && dir.exists(object)) {
    files <- list.files(object, pattern = "\\.Rdata$|\\.rds$",
                        full.names = TRUE, ignore.case = TRUE)
    if (length(files) == 0L)
      stop("no .Rdata/.rds files in ", object)
    m <- regexec(pattern, basename(files), ignore.case = TRUE)
    hit <- regmatches(basename(files), m)
    ok  <- lengths(hit) > 0L
    if (!any(ok))
      stop("no files matching ", pattern, " in ", object,
           "; pass a named character vector of paths instead")
    catch <- vapply(hit[ok], function(h) h[2L], character(1))
    object <- setNames(files[ok], catch)
  }

  if (is.null(names(object)) || any(!nzchar(names(object))))
    stop("'object' must be a named character vector of paths, ",
         "or a directory of *_\\{catch\\}t.Rdata files")

  catch <- suppressWarnings(as.numeric(names(object)))
  if (any(is.na(catch)))
    stop("names(object) must be catch levels coercible to numeric")

  # stable order by catch
  object <- object[order(catch)]

  out <- lapply(object, function(p) {
    if (!file.exists(p) && !dir.exists(p))
      stop("path not found: ", p)
    # ssLoad lives in trajectory.R
    ssLoad(p)
  })
  names(out) <- names(object)
  out
})

#' @rdname getRuns
#' @export
setMethod("getRuns", signature(object = "list"),
  function(object, ...) {

  if (length(object) == 0L) stop("'object' is empty")
  catch <- suppressWarnings(as.numeric(names(object)))
  if (length(catch) == 0L || any(is.na(catch)))
    stop("names(object) must be catch levels coercible to numeric")

  isSS <- function(o)
    is.list(o) && all(c("timeseries", "derived_quants") %in% names(o))
  bad <- !vapply(object, isSS, logical(1))
  if (any(bad))
    stop("list elements must be SS_output-shaped; bad: ",
         paste(names(object)[bad], collapse = ", "))

  object[order(catch)]
})

#' Analytic Kobe II strategy matrix
#'
#' Builds the three probability tables of Liniers et al. (SCRS/2026/178) —
#' \eqn{P(F \\le F_{MSY})}, \eqn{P(SSF \\ge SSF_{MSY})}, and the green quadrant —
#' without Monte Carlo, from Stock Synthesis delta-method moments.
#'
#' @param object A named list of runs (as from \code{\link{getRuns}}) or a
#'   moment grid from \code{\link{kobeMomentGrid}}.
#' @param years Years to tabulate; defaults to the Rec. 21-09 reporting years
#'   when \code{object} is a run list.
#' @param transform Passed to \code{\link{kobeMomentGrid}}.
#' @param asPercent Logical; if \code{TRUE}, return percentages rounded to
#'   integers as in the published table.
#' @param ... Unused.
#' @return A list with elements \code{pNoOverfishing}, \code{pNotOverfished},
#'   \code{green} (each a catch x year matrix) and \code{long} (tidy form).
#' @examples
#' \dontrun{
#' runs <- getRuns("P:/.../SSoutput")
#' k2   <- makeK2SM(runs, years = c(2027, 2030, 2040, 2050, 2070))
#' k2$green
#' }
#' @seealso \code{\link{kobeMomentGrid}}, \code{\link{kobeProbs}},
#'   \code{\link{getRuns}}
#' @export
setGeneric("makeK2SM", function(object, ...) standardGeneric("makeK2SM"))

#' @rdname makeK2SM
#' @export
setMethod("makeK2SM", signature(object = "list"),
  function(object, years = c(2027, 2028, 2029, 2030, 2035, 2040, 2045,
                             2050, 2055, 2060, 2065, 2070),
           transform = c("exact", "current"), asPercent = FALSE, ...) {

  transform <- match.arg(transform)
  grid <- kobeMomentGrid(object, years, transform = transform)
  makeK2SM(grid, years = years, asPercent = asPercent)
})

#' @rdname makeK2SM
#' @export
setMethod("makeK2SM", signature(object = "data.frame"),
  function(object, years = NULL, asPercent = FALSE, ...) {

  need <- c("catch", "year", "muB", "muF", "sigmaB", "sigmaF", "rhoLog")
  if (!all(need %in% names(object)))
    stop("'object' must be kobeMomentGrid output; missing: ",
         paste(setdiff(need, names(object)), collapse = ", "))

  if (is.null(years)) years <- sort(unique(object$year))
  g <- object[object$year %in% years, , drop = FALSE]

  long <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
    r <- g[i, ]
    q <- function(w)
      kobeGreen(r$muB, r$muF, r$sigmaB, r$sigmaF, r$rhoLog, quadrant = w)
    data.frame(catch = r$catch, year = r$year,
               pNoOverfishing = stats::pnorm(0, r$muF, r$sigmaF),
               pNotOverfished = 1 - stats::pnorm(0, r$muB, r$sigmaB),
               green = q("green"),
               stringsAsFactors = FALSE)
  }))
  rownames(long) <- NULL

  pivot <- function(col) {
    m <- reshape(long[, c("catch", "year", col)],
                 idvar = "catch", timevar = "year", direction = "wide")
    names(m) <- sub(paste0("^", col, "\\."), "", names(m))
    rownames(m) <- NULL
    if (asPercent) {
      pct <- m
      pct[-1] <- round(100 * as.matrix(m[-1]))
      return(pct)
    }
    m
  }

  list(pNoOverfishing = pivot("pNoOverfishing"),
       pNotOverfished = pivot("pNotOverfished"),
       green          = pivot("green"),
       long           = long)
})

#' Rebuild reference points from a recovery-time curve
#'
#' Ports the RebuildAge response metrics: rebuilding time (in generation times)
#' from \eqn{B_{lim}} and \eqn{B_{trigger}}, and the biomass fractions that
#' recover in one and half a generation time.
#'
#' @param object A \code{data.frame} with columns \code{initial} (starting
#'   biomass fraction) and \code{TRecover} (years to recover), optionally with
#'   \code{gt}; or an \code{FLBRP} / rebuilding trajectory that
#'   \code{\link{rebuildTime}} understands.
#' @param gt Generation time in years. Used to express times in generation
#'   units; taken from \code{object$gt} when present, otherwise required.
#' @param blim Biomass fraction defining \eqn{B_{lim}} (default 0.30).
#' @param btrigger Biomass fraction defining \eqn{B_{trigger}} (default 0.46).
#' @param ... Passed to \code{\link{rebuildTime}} for \code{FLBRP} methods.
#' @return A one-row \code{data.frame} with \code{Tlim}, \code{Ttrig},
#'   \code{Bgt}, \code{Bgt50} and \code{gt}.
#' @references
#' Kell, L.T. et al. Age-structured rebuilding reference points (RebuildAge).
#' @seealso \code{\link{rebuildTime}}, \code{\link{kobeRebuildCurve}}
#' @export
setGeneric("makeRebuildRefs",
           function(object, ...) standardGeneric("makeRebuildRefs"))

#' @rdname makeRebuildRefs
#' @export
setMethod("makeRebuildRefs", signature(object = "data.frame"),
  function(object, gt = NULL, blim = 0.30, btrigger = 0.46, ...) {

  if (!all(c("initial", "TRecover") %in% names(object)))
    stop("'object' must have columns 'initial' and 'TRecover'")

  if (is.null(gt)) {
    if ("gt" %in% names(object) && any(is.finite(object$gt)))
      gt <- mean(object$gt, na.rm = TRUE)
    else
      stop("supply 'gt' (generation time in years)")
  }
  if (!is.finite(gt) || gt <= 0) stop("'gt' must be a positive finite number")

  ok <- is.finite(object$initial) & is.finite(object$TRecover)
  x  <- object$initial[ok]
  y  <- object$TRecover[ok] / gt
  o  <- order(x); x <- x[o]; y <- y[o]
  keep <- !duplicated(x); x <- x[keep]; y <- y[keep]

  if (length(x) < 2L)
    return(data.frame(Tlim = NA_real_, Ttrig = NA_real_,
                      Bgt = NA_real_, Bgt50 = NA_real_, gt = gt))

  approx_xy <- function(xout)
    as.numeric(stats::approx(x, y, xout = xout, rule = 1)$y)
  approx_yx <- function(yout)
    as.numeric(stats::approx(y, x, xout = yout, rule = 1)$y)

  data.frame(Tlim  = approx_xy(blim),
             Ttrig = approx_xy(btrigger),
             Bgt   = approx_yx(1.00),
             Bgt50 = approx_yx(0.50),
             gt    = gt)
})

#' @rdname makeRebuildRefs
#' @export
setMethod("makeRebuildRefs", signature(object = "FLBRP"),
  function(object, gt = NULL, blim = 0.30, btrigger = 0.46, ...) {

  rb <- rebuild(object, ...)
  rt <- rebuildTime(rb, ...)
  if (!all(c("initial", "TRecover") %in% names(rt))) {
    # rebuildTime naming varies; normalise
    nms <- names(rt)
    if (all(c("initial", "year") %in% nms))
      names(rt)[match("year", nms)] <- "TRecover"
    else if (ncol(rt) >= 2L)
      names(rt)[1:2] <- c("initial", "TRecover")
  }
  if (is.null(gt)) {
    gt <- tryCatch(c(genTimeFn(object)), error = function(e) NA_real_)
    if (!is.finite(gt))
      stop("could not compute generation time; pass gt= explicitly")
  }
  makeRebuildRefs(as.data.frame(rt), gt = gt, blim = blim, btrigger = btrigger)
})
