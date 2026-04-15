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
