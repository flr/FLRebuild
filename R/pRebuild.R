# ------------------------------------------------------------------------------
# pRebuild : rebuilding diagnostics from an FLBRP object
#
# Returns an FLPar (params x iter) with rows:
#   B0, Bmsy, Fmsy, Fcrash, GT, SPR0,
#   SSB_pR, SSB_pR_frac,
#   Trec_F0, Trec_Fmsy, Trec_HCR,
#   Trec_F0_GT, Trec_Fmsy_GT, Trec_HCR_GT,
#   Bstar_F0, Bstar_Fmsy,
#   Bstar_F0_frac, Bstar_Fmsy_frac
#
# Projection engine: age-structured Baranov catch equation with recruitment
# from the FLBRP's stock-recruit function. No FLasher::fwd, no FLRebuild -- so
# it is version-independent.
#
# Methods:
#   pRebuild(FLBRP)                      -- single stock, iter-aware
#   pRebuild(FLBRPs, parallel = TRUE)    -- many stocks, optional parallel
# ------------------------------------------------------------------------------

library(FLCore); library(FLBRP)

# ==============================================================================
# GENERICS AND METHODS
# ==============================================================================

setGeneric("pRebuild", function(object, ...) standardGeneric("pRebuild"))

## ---- FLBRP method: single stock, loops over iters ---------------------------
setMethod("pRebuild", signature(object = "FLBRP"),
function(object, pRec = 0.20, maxYears = 200) {

  ni <- dims(object)$iter
  param_names <- c("B0","Bmsy","Fmsy","Fcrash","GT","SPR0",
                   "SSB_pR","SSB_pR_frac",
                   "Trec_F0","Trec_Fmsy","Trec_HCR",
                   "Trec_F0_GT","Trec_Fmsy_GT","Trec_HCR_GT",
                   "Bstar_F0","Bstar_Fmsy",
                   "Bstar_F0_frac","Bstar_Fmsy_frac")

  mat <- matrix(NA_real_, nrow = length(param_names), ncol = ni,
                dimnames = list(params = param_names, iter = seq_len(ni)))

  for (i in seq_len(ni)) {
    obj_i <- iter(object, i)
    mat[, i] <- tryCatch(
      .pRebuild_one(obj_i, pRec = pRec, maxYears = maxYears),
      error = function(e) {
        message("pRebuild iter ", i, " failed: ", conditionMessage(e))
        rep(NA_real_, length(param_names))
      })
  }

  FLPar(mat, units = "NA")
})

## ---- FLBRPs method: many stocks, optional parallel --------------------------
setMethod("pRebuild", signature(object = "FLBRPs"),
function(object,
         pRec     = 0.20,
         maxYears = 200,
         parallel = FALSE,
         ncores   = NULL,
         combine  = c("list", "flpar")) {

  combine <- match.arg(combine)
  stocks  <- names(object)
  if (is.null(stocks)) stocks <- paste0("stk", seq_along(object))

  if (parallel) {
    if (!requireNamespace("foreach",    quietly = TRUE) ||
        !requireNamespace("doParallel", quietly = TRUE))
      stop("Install foreach and doParallel for parallel = TRUE")

    if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)
    ncores <- min(ncores, length(object))

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    `%dopar%` <- foreach::`%dopar%`
    res <- foreach::foreach(i = seq_along(object),
                            .packages = c("FLCore","FLBRP","methods"),
                            .errorhandling = "pass") %dopar% {
      pRebuild(object[[i]], pRec = pRec, maxYears = maxYears)
    }
  } else {
    res <- lapply(seq_along(object), function(i)
      pRebuild(object[[i]], pRec = pRec, maxYears = maxYears))
  }

  names(res) <- stocks

  if (combine == "list") return(res)

  ## stack single-iter FLPars column-wise into one FLPar (cols = stocks)
  ok  <- vapply(res, function(x) is(x, "FLPar"), logical(1))
  if (!any(ok)) return(res)
  mat <- do.call(cbind, lapply(res[ok], function(x) x@.Data))
  colnames(mat) <- stocks[ok]
  FLPar(mat, units = "NA")
})


# ==============================================================================
# CORE SINGLE-ITER ENGINE
# ==============================================================================

.pRebuild_one <- function(brp, pRec, maxYears) {

  rp     <- refpts(brp)
  B0     <- c(rp["virgin", "ssb"])
  Bmsy   <- c(rp["msy",    "ssb"])
  Fmsy   <- c(rp["msy",    "harvest"])
  Fcrash <- c(rp["crash",  "harvest"])
  GT     <- c(gt(brp))
  SPR0   <- tryCatch(c(spr0(brp)), error = function(e) NA_real_)

  SSB_pR <- .ssb_at_p_recruit(brp, pRec, B0)
  init   <- pRec * B0

  T_F0   <- .safe_time(brp, init, Bmsy, targetF = 0,     maxYears)
  T_Fmsy <- .safe_time(brp, init, Bmsy, targetF = Fmsy,  maxYears)
  T_HCR  <- .safe_time(brp, init, Bmsy, targetF = "hcr", maxYears, Fmsy = Fmsy)

  Bstar_F0   <- .bstar_one_gen(brp, GT, Bmsy, B0, targetF = 0,    Fmsy = Fmsy)
  Bstar_Fmsy <- .bstar_one_gen(brp, GT, Bmsy, B0, targetF = Fmsy, Fmsy = Fmsy)

  c(B0              = B0,
    Bmsy            = Bmsy,
    Fmsy            = Fmsy,
    Fcrash          = Fcrash,
    GT              = GT,
    SPR0            = SPR0,
    SSB_pR          = SSB_pR,
    SSB_pR_frac     = SSB_pR / B0,
    Trec_F0         = T_F0,
    Trec_Fmsy       = T_Fmsy,
    Trec_HCR        = T_HCR,
    Trec_F0_GT      = T_F0   / GT,
    Trec_Fmsy_GT    = T_Fmsy / GT,
    Trec_HCR_GT     = T_HCR  / GT,
    Bstar_F0        = Bstar_F0,
    Bstar_Fmsy      = Bstar_Fmsy,
    Bstar_F0_frac   = Bstar_F0   / B0,
    Bstar_Fmsy_frac = Bstar_Fmsy / B0)
}


# ==============================================================================
# HELPERS
# ==============================================================================

## ---- SSB at which recruitment is p * virgin ---------------------------------
.ssb_at_p_recruit <- function(brp, p, B0) {
  nm <- SRModelName(model(brp))
  if (grepl("bevholt", nm, ignore.case = TRUE)) {
    b <- c(params(brp)["b"]); return(p * b / (1 - p))
  }
  if (grepl("segreg|hockey", nm, ignore.case = TRUE)) {
    b <- c(params(brp)["b"]); return(p * b)
  }
  sr_pred <- .make_sr(brp)
  target  <- p * sr_pred(B0)
  tryCatch(uniroot(function(s) sr_pred(s) - target, c(1e-6, B0))$root,
           error = function(e) NA_real_)
}


## ---- error-tolerant wrapper around .rebuild_time ----------------------------
.safe_time <- function(...) {
  tryCatch(.rebuild_time(...),
           error = function(e) {
             message("rebuild_time failed: ", conditionMessage(e))
             NA_real_
           })
}


## ---- SR predictor closure ---------------------------------------------------
.make_sr <- function(brp) {
  nm <- SRModelName(model(brp))
  a  <- tryCatch(c(params(brp)["a"]), error = function(e) NA_real_)
  b  <- tryCatch(c(params(brp)["b"]), error = function(e) NA_real_)

  if (!is.na(a) && !is.na(b)) {
    if (grepl("bevholt", nm, ignore.case = TRUE))
      return(function(ssb) a * ssb / (b + ssb))
    if (grepl("ricker",  nm, ignore.case = TRUE))
      return(function(ssb) a * ssb * exp(-b * ssb))
    if (grepl("segreg|hockey", nm, ignore.case = TRUE))
      return(function(ssb) ifelse(ssb >= b, a * b, a * ssb))
  }
  ## last-resort fallback via predictModel
  function(ssb) {
    v <- tryCatch(c(predict(as(brp, "FLSR"), ssb = FLQuant(ssb))),
                  error = function(e) NA_real_)
    if (is.na(v)) 0 else v
  }
}


## ---- Baranov projection: recovery-time engine -------------------------------
.rebuild_time <- function(brp, init, target, targetF, maxYears, Fmsy = NULL) {

  eql <- brp(brp)

  M   <- c(m(eql)[, 1])
  mat <- c(mat(eql)[, 1])
  swt <- c(stock.wt(eql)[, 1])

  ## selectivity waterfall: harvest -> catch.sel -> flat
  sel <- c(harvest(eql)[, 1])
  if (all(is.na(sel)) || max(sel, na.rm = TRUE) <= 0)
    sel <- tryCatch(c(catch.sel(brp)[, 1]), error = function(e) rep(NA, length(M)))
  if (all(is.na(sel)) || max(sel, na.rm = TRUE) <= 0)
    sel <- rep(1, length(M))
  else
    sel <- sel / max(sel, na.rm = TRUE)

  nage    <- length(M)
  sr_pred <- .make_sr(brp)

  ## initial N-at-age scaled so year-1 SSB = init
  Neq    <- c(stock.n(eql)[, 1])
  ssb_eq <- c(ssb(eql)[, 1])
  if (is.na(ssb_eq) || ssb_eq <= 0) return(NA_real_)
  N       <- Neq * (init / ssb_eq)
  ssb_now <- sum(N * mat * swt, na.rm = TRUE)
  if (ssb_now >= target) return(0L)         # already recovered

  for (y in 1:maxYears) {
    Fy <- if (identical(targetF, "hcr"))
            Fmsy * min(1, ssb_now / target)
          else
            as.numeric(targetF)

    Fa <- sel * Fy
    Za <- Fa + M

    R  <- sr_pred(ssb_now)
    if (is.na(R) || R < 0) R <- 0

    Nn <- numeric(nage)
    Nn[1]      <- R
    if (nage >= 2)
      Nn[2:nage] <- N[1:(nage - 1)] * exp(-Za[1:(nage - 1)])
    Nn[nage]   <- Nn[nage] + N[nage] * exp(-Za[nage])   # plus group

    N       <- Nn
    ssb_now <- sum(N * mat * swt, na.rm = TRUE)

    if (is.na(ssb_now)) return(NA_real_)
    if (ssb_now >= target) return(y)
  }
  NA_real_
}


## ---- bisection: starting SSB from which recovery in 1 GT is possible --------
.bstar_one_gen <- function(brp, GT, Bmsy, B0, targetF, Fmsy) {

  horizon <- max(50, ceiling(3 * GT))

  f <- function(frac) {
    tt <- .safe_time(brp, init = frac * B0, target = Bmsy,
                     targetF = targetF, maxYears = horizon,
                     Fmsy = Fmsy)
    if (is.na(tt)) return(NA_real_)
    tt - GT                          # positive if too slow; negative if fast enough
  }

  lo <- 0.02
  hi <- (Bmsy / B0) * 0.999

  fl <- f(lo)
  fh <- f(hi)

  if (is.na(fl) || is.na(fh)) return(NA_real_)

  ## normal case: fl > 0 (slow from lo), fh < 0 (fast from hi) -> bisect
  if (fl > 0 && fh < 0) {
    root <- tryCatch(
      uniroot(f, c(lo, hi), tol = 1e-3)$root,
      error = function(e) NA_real_)
    return(if (is.na(root)) NA_real_ else root * B0)
  }
  ## degenerate cases
  if (fl <= 0) return(lo * B0)       # already fast enough from lo
  if (fh >  0) return(hi * B0)       # cannot make 1GT even from ~Bmsy
  NA_real_
}
