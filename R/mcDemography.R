#' Monte Carlo Demographic Estimates
#'
#' Generic interface to compute deterministic and stochastic demographic
#' quantities from stock objects.
#'
#' @param object A stock object.
#' @param ... Additional method-specific arguments.
#'
#' @return A list with deterministic estimates, stochastic draws, and summary.
#' @export
setGeneric("mcDemography", function(object, ...) standardGeneric("mcDemography"))

.demographyFromLeslie <- function(L) {
  eig <- eigen(L)
  vals <- eig$values
  mods <- sort(Mod(vals), decreasing=TRUE)
  damp <- if (length(mods) > 1L) mods[1] / mods[2] else NA_real_

  lambda <- Re(vals[which.max(Mod(vals))])
  r <- log(lambda)

  fec <- L[1, ]
  ages <- length(fec)
  lx <- rep(0, ages)
  lx[1] <- 1
  if (ages > 1) {
    for (i in 2:ages)
      lx[i] <- lx[i - 1] * L[i, i - 1]
  }

  R0 <- sum(lx * fec, na.rm=TRUE)
  gt <- if (is.finite(R0) && R0 > 0) sum((1:ages) * lx * fec, na.rm=TRUE) / R0 else NA_real_

  c(lambda=lambda, r=r, R0=R0, G=gt, damp=damp)
}

.buildLeslieFlStock <- function(object,
                                fbar=0,
                                femalePupsPerYear=1.0) {
  ageNames <- dimnames(FLCore::stock.n(object))$age
  ages <- as.numeric(ageNames)
  nages <- length(ages)

  mAtAge <- as.numeric(FLCore::m(object)[, 1, 1, 1, 1, 1])
  if (length(mAtAge) != nages)
    stop("Could not extract age-specific natural mortality from FLStock")

  matAtAge <- as.numeric(FLCore::mat(object)[, 1, 1, 1, 1, 1])
  fec <- pmax(matAtAge, 0) * femalePupsPerYear

  surv <- exp(-(mAtAge + fbar))
  if (length(surv) > 1L)
    surv <- surv[1:(length(surv) - 1L)]

  L <- matrix(0, nrow=nages, ncol=nages)
  L[1, ] <- fec
  if (nages > 1)
    L[cbind(2:nages, 1:(nages - 1))] <- surv

  L
}

.summariseDraws <- function(drawsDf, probs=c(0.025, 0.5, 0.975)) {
  out <- lapply(drawsDf, function(x) {
    stats::quantile(x, probs=probs, na.rm=TRUE, names=FALSE)
  })
  q <- do.call(rbind, out)
  colnames(q) <- paste0("q", probs * 100)
  data.frame(metric=rownames(q), q, row.names=NULL, check.names=FALSE)
}

#' @rdname mcDemography
#' @param fbar Fishing mortality used in the Leslie matrix.
#' @param n Number of stochastic simulations.
#' @param mCv Lognormal coefficient of variation applied to mortality.
#' @param fecCv Lognormal coefficient of variation applied to fecundity.
#' @param seed Optional RNG seed for reproducibility.
#' @exportMethod mcDemography
setMethod("mcDemography", signature(object="FLBRP"),
function(object,
         fbar=FLCore::FLQuant(c(FLCore::refpts(object)["msy", "harvest"])),
         n=2000L,
         mCv=0.2,
         fecCv=0.2,
         seed=NULL,
         ...) {

  if (!is.null(seed))
    set.seed(seed)

  L <- FLife::leslie(object, fbar=fbar)[drop=TRUE]
  det <- .demographyFromLeslie(L)

  n <- as.integer(n)
  draws <- matrix(NA_real_, nrow=n, ncol=length(det))
  colnames(draws) <- names(det)

  for (i in seq_len(n)) {
    mMult <- rlnorm(1, meanlog=-0.5 * log(1 + mCv^2), sdlog=sqrt(log(1 + mCv^2)))
    fMult <- rlnorm(1, meanlog=-0.5 * log(1 + fecCv^2), sdlog=sqrt(log(1 + fecCv^2)))
    Ls <- L
    if (nrow(Ls) > 1)
      Ls[2:nrow(Ls), ] <- Ls[2:nrow(Ls), ]^mMult
    Ls[1, ] <- Ls[1, ] * fMult
    draws[i, ] <- .demographyFromLeslie(Ls)
  }

  drawsDf <- as.data.frame(draws)
  list(
    deterministic=det,
    stochastic=drawsDf,
    summary=.summariseDraws(drawsDf)
  )
})

#' @rdname mcDemography
#' @param femalePupsPerYear Scalar fecundity scale for mature ages.
#' @exportMethod mcDemography
setMethod("mcDemography", signature(object="FLStock"),
function(object,
         fbar=0,
         femalePupsPerYear=1.0,
         n=2000L,
         mCv=0.2,
         fecCv=0.2,
         seed=NULL,
         ...) {

  if (!is.null(seed))
    set.seed(seed)

  L <- .buildLeslieFlStock(
    object=object,
    fbar=fbar,
    femalePupsPerYear=femalePupsPerYear
  )
  det <- .demographyFromLeslie(L)

  n <- as.integer(n)
  draws <- matrix(NA_real_, nrow=n, ncol=length(det))
  colnames(draws) <- names(det)

  for (i in seq_len(n)) {
    mMult <- rlnorm(1, meanlog=-0.5 * log(1 + mCv^2), sdlog=sqrt(log(1 + mCv^2)))
    fMult <- rlnorm(1, meanlog=-0.5 * log(1 + fecCv^2), sdlog=sqrt(log(1 + fecCv^2)))
    Ls <- L
    if (nrow(Ls) > 1)
      Ls[2:nrow(Ls), ] <- Ls[2:nrow(Ls), ]^mMult
    Ls[1, ] <- Ls[1, ] * fMult
    draws[i, ] <- .demographyFromLeslie(Ls)
  }

  drawsDf <- as.data.frame(draws)
  list(
    deterministic=det,
    stochastic=drawsDf,
    summary=.summariseDraws(drawsDf)
  )
})
