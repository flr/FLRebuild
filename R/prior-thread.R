#' Build Canonical Prior Table
#'
#' Build a canonical prior table (`param`, `mean`, `sd`, `cv`, `source`) from
#' either an `FLStock` plus `jabbaPriors()` method, or an existing JABBA prior
#' object (`list(pr, pr.sd)`).
#'
#' @param object An `FLStock` or a JABBA prior object with `pr` and `pr.sd`.
#' @param method Prior pathway used when `object` is `FLStock`.
#' @param source Label stored in canonical output.
#' @param priorCv Default coefficient of variation used by `jabbaPriors()`.
#' @param ... Additional arguments forwarded to `jabbaPriors()` when relevant.
#'
#' @return A data.frame with columns `param`, `mean`, `sd`, `cv`, `source`.
#' @export
canonicalPriors <- function(object,
                            method=c("ices", "fishlife", "flbrp"),
                            source=NULL,
                            priorCv=0.30,
                            ...) {
  method <- match.arg(method)

  priors <- NULL
  if (methods::is(object, "FLStock")) {
    priors <- jabbaPriors(object=object, method=method, prior.cv=priorCv, ...)
    if (is.null(source))
      source <- method
  } else if (is.list(object) && all(c("pr", "pr.sd") %in% names(object))) {
    priors <- object
    if (is.null(source))
      source <- "provided"
  } else {
    stop("'object' must be an FLStock or a prior list with 'pr' and 'pr.sd'")
  }

  .validateJabbaPriors(priors)

  param <- dimnames(priors$pr)$params
  meanVals <- as.numeric(c(priors$pr))
  sdVals <- as.numeric(c(priors$pr.sd))
  names(meanVals) <- param
  names(sdVals) <- dimnames(priors$pr.sd)$params
  sdVals <- sdVals[match(param, names(sdVals))]
  cvVals <- ifelse(is.finite(meanVals) & meanVals != 0, sdVals / meanVals, NA_real_)

  data.frame(
    param=param,
    mean=as.numeric(meanVals),
    sd=as.numeric(sdVals),
    cv=as.numeric(cvVals),
    source=rep(as.character(source), length(param)),
    stringsAsFactors=FALSE
  )
}

#' Build Prior Covariance Matrix
#'
#' Build a covariance matrix from canonical prior SDs and an optional
#' correlation matrix.
#'
#' @param canonical Canonical prior data.frame from `canonicalPriors()`.
#' @param correlation Optional correlation matrix. If `NULL`, independence is
#'   assumed (diagonal covariance).
#'
#' @return Named covariance matrix.
#' @export
priorCovariance <- function(canonical, correlation=NULL) {
  req <- c("param", "sd")
  if (!is.data.frame(canonical) || !all(req %in% names(canonical)))
    stop("'canonical' must be a data.frame containing: ", paste(req, collapse=", "))

  pars <- as.character(canonical$param)
  sds <- as.numeric(canonical$sd)
  if (any(!is.finite(sds)))
    stop("'canonical$sd' must be finite values")

  sdDiag <- diag(sds)
  dimnames(sdDiag) <- list(pars, pars)

  if (is.null(correlation)) {
    covMat <- sdDiag %*% diag(length(sds)) %*% sdDiag
    return(covMat)
  }

  if (!is.matrix(correlation))
    stop("'correlation' must be a matrix")
  if (nrow(correlation) != length(pars) || ncol(correlation) != length(pars))
    stop("'correlation' dimensions must match number of parameters in canonical")

  if (is.null(rownames(correlation)) || is.null(colnames(correlation))) {
    rownames(correlation) <- pars
    colnames(correlation) <- pars
  }

  correlation <- correlation[pars, pars, drop=FALSE]
  if (any(abs(diag(correlation) - 1) > 1e-8))
    stop("Diagonal of 'correlation' must be 1")
  if (max(abs(correlation - t(correlation)), na.rm=TRUE) > 1e-8)
    stop("'correlation' must be symmetric")

  covMat <- sdDiag %*% correlation %*% sdDiag
  dimnames(covMat) <- list(pars, pars)
  covMat
}

#' Convert Canonical Priors to JABBA Prior Object
#'
#' @param canonical Canonical prior data.frame.
#'
#' @return A JABBA prior object with `pr` and `pr.sd` as `FLPar`.
#' @export
asJabbaPriors <- function(canonical) {
  req <- c("param", "mean", "sd")
  if (!is.data.frame(canonical) || !all(req %in% names(canonical)))
    stop("'canonical' must be a data.frame containing: ", paste(req, collapse=", "))

  meanVals <- as.numeric(canonical$mean)
  sdVals <- as.numeric(canonical$sd)
  names(meanVals) <- canonical$param
  names(sdVals) <- canonical$param

  priors <- list(
    pr=FLCore::FLPar(meanVals),
    pr.sd=FLCore::FLPar(sdVals)
  )
  .validateJabbaPriors(priors)
  priors
}

#' Convert Canonical Priors to SPiCT Input Template
#'
#' @param canonical Canonical prior data.frame.
#' @param covariance Optional covariance matrix aligned to `canonical$param`.
#'
#' @return A named list template for SPiCT prior mapping.
#' @export
asSpictPriors <- function(canonical, covariance=NULL) {
  req <- c("param", "mean", "sd", "cv")
  if (!is.data.frame(canonical) || !all(req %in% names(canonical)))
    stop("'canonical' must be a data.frame containing: ", paste(req, collapse=", "))

  meanVals <- as.numeric(canonical$mean)
  sdVals <- as.numeric(canonical$sd)
  cvVals <- as.numeric(canonical$cv)
  names(meanVals) <- canonical$param
  names(sdVals) <- canonical$param
  names(cvVals) <- canonical$param

  list(
    mean=meanVals,
    sd=sdVals,
    cv=cvVals,
    covariance=covariance,
    notes="Map these canonical values to SPiCT parameterization."
  )
}

#' Convert Canonical Priors to FLICC Input Template
#'
#' @param canonical Canonical prior data.frame.
#' @param covariance Optional covariance matrix aligned to `canonical$param`.
#'
#' @return A named list template for FLICC prior mapping.
#' @export
asFliccPriors <- function(canonical, covariance=NULL) {
  req <- c("param", "mean", "sd", "cv")
  if (!is.data.frame(canonical) || !all(req %in% names(canonical)))
    stop("'canonical' must be a data.frame containing: ", paste(req, collapse=", "))

  meanVals <- as.numeric(canonical$mean)
  sdVals <- as.numeric(canonical$sd)
  cvVals <- as.numeric(canonical$cv)
  names(meanVals) <- canonical$param
  names(sdVals) <- canonical$param
  names(cvVals) <- canonical$param

  list(
    mean=meanVals,
    sd=sdVals,
    cv=cvVals,
    covariance=covariance,
    notes="Map these canonical values to FLICC input assumptions."
  )
}

#' Compare Prior Approaches for a Stock
#'
#' Build and compare canonical priors across multiple `jabbaPriors()` methods
#' in a single call.
#'
#' @param object An `FLStock`.
#' @param methods Character vector of methods to compare.
#' @param priorCv Default coefficient of variation passed to `jabbaPriors()`.
#' @param strict If `TRUE`, stop on the first method error. If `FALSE`, keep
#'   successful methods and attach failures in `failures`.
#' @param ... Additional arguments forwarded to `jabbaPriors()`.
#'
#' @return A list with:
#' \describe{
#'   \item{canonicalByMethod}{Named list of canonical prior data.frames.}
#'   \item{comparison}{Single combined comparison table.}
#'   \item{failures}{Named list of error messages for failed methods (if any).}
#' }
#' @export
comparePriorApproaches <- function(object,
                                   methods=c("ices", "fishlife", "flbrp"),
                                   priorCv=0.30,
                                   strict=TRUE,
                                   ...) {
  if (!methods::is(object, "FLStock"))
    stop("'object' must be an FLStock")
  methods <- as.character(methods)
  if (length(methods) == 0L)
    stop("'methods' must contain at least one method")

  canonicalByMethod <- list()
  failures <- list()

  for (mth in methods) {
    out <- try(
      canonicalPriors(
        object=object,
        method=mth,
        source=mth,
        priorCv=priorCv,
        ...
      ),
      silent=TRUE
    )
    if (inherits(out, "try-error")) {
      msg <- conditionMessage(attr(out, "condition"))
      failures[[mth]] <- msg
      if (isTRUE(strict))
        stop("Method '", mth, "' failed: ", msg)
    } else {
      canonicalByMethod[[mth]] <- out
    }
  }

  if (length(canonicalByMethod) == 0L)
    stop("No methods succeeded")

  comparison <- do.call(rbind, canonicalByMethod)
  rownames(comparison) <- NULL

  list(
    canonicalByMethod=canonicalByMethod,
    comparison=comparison,
    failures=failures
  )
}

.extractAgeVector <- function(x) {
  vals <- as.numeric(x[, 1, 1, 1, 1, 1])
  vals[is.na(vals)] <- 0
  vals
}

.yprCurveFromStock <- function(object,
                               fRange=c(0, 2),
                               nF=101L,
                               selectivity=NULL) {
  mAtAge <- .extractAgeVector(FLCore::m(object))
  wtAtAge <- .extractAgeVector(FLCore::stock.wt(object))
  matAtAge <- .extractAgeVector(FLCore::mat(object))
  mSpwn <- .extractAgeVector(FLCore::m.spwn(object))
  hSpwn <- .extractAgeVector(FLCore::harvest.spwn(object))

  nAges <- min(length(mAtAge), length(wtAtAge), length(matAtAge), length(mSpwn), length(hSpwn))
  mAtAge <- mAtAge[seq_len(nAges)]
  wtAtAge <- wtAtAge[seq_len(nAges)]
  matAtAge <- matAtAge[seq_len(nAges)]
  mSpwn <- mSpwn[seq_len(nAges)]
  hSpwn <- hSpwn[seq_len(nAges)]

  if (is.null(selectivity)) {
    sel <- .extractAgeVector(FLCore::harvest(object))
    sel <- sel[seq_len(nAges)]
    if (!any(is.finite(sel)) || max(sel, na.rm=TRUE) <= 0)
      sel <- rep(1, nAges)
    sel <- pmax(sel, 0)
    sel <- sel / max(sel, na.rm=TRUE)
  } else {
    sel <- as.numeric(selectivity)
    if (length(sel) != nAges)
      stop("length(selectivity) must match number of age classes")
    if (max(sel, na.rm=TRUE) <= 0)
      stop("selectivity must contain at least one positive value")
    sel <- pmax(sel, 0) / max(sel, na.rm=TRUE)
  }

  fGrid <- seq(fRange[1], fRange[2], length.out=as.integer(nF))
  ypr <- rep(NA_real_, length(fGrid))
  spr <- rep(NA_real_, length(fGrid))

  for (i in seq_along(fGrid)) {
    fVal <- fGrid[i]
    fAtAge <- fVal * sel
    zAtAge <- mAtAge + fAtAge

    lAtAge <- rep(NA_real_, nAges)
    lAtAge[1] <- 1
    if (nAges > 1) {
      for (a in 2:nAges)
        lAtAge[a] <- lAtAge[a - 1] * exp(-zAtAge[a - 1])
    }

    cAtAge <- lAtAge * (fAtAge / pmax(zAtAge, 1e-12)) * (1 - exp(-zAtAge))
    ypr[i] <- sum(cAtAge * wtAtAge, na.rm=TRUE)

    sprSpawn <- lAtAge * exp(-(mAtAge * mSpwn + fAtAge * hSpwn))
    spr[i] <- sum(sprSpawn * matAtAge * wtAtAge, na.rm=TRUE)
  }

  data.frame(
    f=fGrid,
    ypr=ypr,
    spr=spr,
    stringsAsFactors=FALSE
  )
}

#' Build FLICC Priors from Yield-Per-Recruit
#'
#' Build a FLICC-oriented prior object from an `FLStock` using a
#' yield-per-recruit (YPR) and spawning-per-recruit (SPR) curve.
#'
#' @param object An `FLStock`.
#' @param fRange Range of fishing mortality used for YPR curve.
#' @param nF Number of F grid points.
#' @param selectivity Optional age selectivity vector. If `NULL`, harvest-at-age
#'   is used and normalized.
#' @param priorCv Coefficient of variation used when deriving SDs for scalar
#'   reference quantities.
#'
#' @return A list with:
#' \describe{
#'   \item{curve}{Data.frame with `f`, `ypr`, and `spr`.}
#'   \item{reference}{Named numeric vector with `fMsyProxy`, `yprMsyProxy`,
#'   `spr0`, `sprAtFmsy`, and `sprRatioAtFmsy`.}
#'   \item{canonical}{Canonical-style data.frame for YPR-derived scalar priors.}
#'   \item{covariance}{Covariance matrix from canonical SDs.}
#' }
#' @export
asFliccPriorsYpr <- function(object,
                             fRange=c(0, 2),
                             nF=101L,
                             selectivity=NULL,
                             priorCv=0.30) {
  if (!methods::is(object, "FLStock"))
    stop("'object' must be an FLStock")
  if (length(fRange) != 2L || fRange[1] < 0 || fRange[2] <= fRange[1])
    stop("'fRange' must be length-2 with 0 <= fRange[1] < fRange[2]")
  nF <- as.integer(nF)
  if (nF < 5L)
    stop("'nF' must be at least 5")

  curve <- .yprCurveFromStock(
    object=object,
    fRange=fRange,
    nF=nF,
    selectivity=selectivity
  )

  iMax <- which.max(curve$ypr)
  spr0 <- curve$spr[1]
  sprAtFmsy <- curve$spr[iMax]
  sprRatio <- sprAtFmsy / spr0

  reference <- c(
    fMsyProxy=curve$f[iMax],
    yprMsyProxy=curve$ypr[iMax],
    spr0=spr0,
    sprAtFmsy=sprAtFmsy,
    sprRatioAtFmsy=sprRatio
  )

  canonical <- data.frame(
    param=names(reference),
    mean=as.numeric(reference),
    sd=as.numeric(reference) * priorCv,
    cv=rep(priorCv, length(reference)),
    source=rep("ypr", length(reference)),
    stringsAsFactors=FALSE
  )
  canonical$sd[!is.finite(canonical$sd)] <- NA_real_
  canonical$sd[canonical$sd < 0] <- abs(canonical$sd[canonical$sd < 0])

  canonicalCov <- canonical[is.finite(canonical$sd), , drop=FALSE]
  covariance <- if (nrow(canonicalCov) > 0L) priorCovariance(canonicalCov) else matrix(numeric(0), 0, 0)

  list(
    curve=curve,
    reference=reference,
    canonical=canonical,
    covariance=covariance
  )
}
