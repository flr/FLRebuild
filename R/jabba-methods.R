#' @rdname jabbaInput
setGeneric("jabbaInput", function(object, ...) standardGeneric("jabbaInput"))

#' @rdname jabbaPriors
setGeneric("jabbaPriors",
  function(object, method=c("ices", "fishlife", "flbrp"), ...)
    standardGeneric("jabbaPriors"))

#' @rdname runJABBA
setGeneric("runJABBA", function(object, ...) standardGeneric("runJABBA"))

.validateJabbaPriors <- function(priors) {
  if (!is.list(priors) || !all(c("pr", "pr.sd") %in% names(priors)))
    stop("priors must be a list with elements 'pr' and 'pr.sd'")

  if (!methods::is(priors$pr, "FLPar"))
    stop("priors$pr must be an FLPar")

  if (!methods::is(priors$pr.sd, "FLPar"))
    stop("priors$pr.sd must be an FLPar")

  req <- c("r", "psi", "shape")
  pnm <- dimnames(priors$pr)$params
  snm <- dimnames(priors$pr.sd)$params

  if (!all(req %in% pnm))
    stop("priors$pr must contain parameters: ", paste(req, collapse=", "))
  if (!all(req %in% snm))
    stop("priors$pr.sd must contain parameters: ", paste(req, collapse=", "))

  opt <- intersect(c("k", "current"), pnm)
  if (length(opt) > 0L && !all(opt %in% snm))
    stop("priors$pr.sd must also contain optional parameters found in priors$pr: ",
         paste(opt, collapse=", "))

  invisible(priors)
}

.validateJabbaInput <- function(catch, index=NULL) {
  if (!is.data.frame(catch))
    stop("'catch' must be a data.frame")
  if (!all(c("year", "catch") %in% names(catch)))
    stop("'catch' must contain columns 'year' and 'catch'")

  if (!is.null(index)) {
    if (!is.data.frame(index))
      stop("'index' must be a data.frame")
    if (!all(c("year", "index") %in% names(index)))
      stop("'index' must contain columns 'year' and 'index'")
  }

  invisible(list(catch=catch, index=index))
}

.stockDepletion <- function(object, b0, initial.yrs=5, current.yrs=5) {
  ts <- tseries(object)
  yrs0 <- seq(min(ts$year), length.out=initial.yrs)
  yrs1 <- seq(max(ts$year) - current.yrs + 1, max(ts$year))
  c(
    psi     = mean(ts$ssb[ts$year %in% yrs0], na.rm=TRUE) / b0,
    current = mean(ts$ssb[ts$year %in% yrs1], na.rm=TRUE) / b0
  )
}

.index2df <- function(x, value.name="index") {
  if (methods::is(x, "FLQuant")) {
    yrs  <- as.integer(dimnames(x)$year)
    vals <- as.numeric(c(x))
    return(data.frame(year=yrs, index=vals))
  }
  if (is.data.frame(x)) {
    if (!all(c("year", value.name) %in% names(x)))
      stop("index data.frame must contain columns 'year' and '", value.name, "'")
    return(x[, c("year", value.name), drop=FALSE])
  }
  stop("Unsupported index object passed to .index2df()")
}

.asIndexDf <- function(x, value.name="index") {
  if (exists(".index2df", mode="function", inherits=TRUE))
    return(.index2df(x, value.name=value.name))

  if (methods::is(x, "FLQuant")) {
    yrs  <- as.integer(dimnames(x)$year)
    vals <- as.numeric(c(x))
    return(data.frame(year=yrs, index=vals))
  }
  if (is.data.frame(x)) {
    if (!all(c("year", value.name) %in% names(x)))
      stop("index data.frame must contain columns 'year' and '", value.name, "'")
    return(x[, c("year", value.name), drop=FALSE])
  }
  stop("Unsupported index object passed to jabbaInput(index=...)")
}

.jabbaPellatParams <- function(fmsy, bmsy, b0, interval=NULL) {
  fmsy <- as.numeric(fmsy)[1]
  bmsy <- as.numeric(bmsy)[1]
  b0   <- as.numeric(b0)[1]
  shape <- bmsy / b0

  if (is.null(interval))
    interval <- if (shape < 0.37) c(0.10, 0.99) else c(1.01, 10.0)

  m <- stats::optimize(
    function(x) abs(shape - (1 / x)^(1 / (x - 1))),
    interval=interval
  )$minimum
  r <- fmsy * (m - 1) / (1 - 1 / m)
  FLCore::FLPar(c(r=r, shape=shape, k=b0))
}

.jabbaPriorObject <- function(pars, sds=NULL, prior.cv=0.30) {
  if (is.list(pars))
    pars <- unlist(pars, use.names=TRUE)
  pars <- stats::setNames(as.numeric(pars), names(pars))

  req <- c("r", "psi", "shape")
  if (!all(req %in% names(pars)))
    stop(
      "Prior means must contain: ", paste(req, collapse=", "),
      ". Found: ", paste(names(pars), collapse=", ")
    )

  keep <- intersect(c("r", "psi", "shape", "k", "current"), names(pars))
  pr <- FLCore::FLPar(pars[keep])

  if (is.null(sds)) {
    pr.sd <- pr %/% pr * prior.cv
  } else {
    if (is.list(sds))
      sds <- unlist(sds, use.names=TRUE)
    sds <- stats::setNames(as.numeric(sds), names(sds))
    if (!all(keep %in% names(sds)))
      stop(
        "Prior sds must contain: ", paste(keep, collapse=", "),
        ". Found: ", paste(names(sds), collapse=", ")
      )
    pr.sd <- FLCore::FLPar(sds[keep])
  }

  list(pr=pr, pr.sd=pr.sd)
}

.named_num <- function(x, nm) {
  stats::setNames(as.numeric(x)[1], nm)
}

.jabbaPriorsICES <- function(object, prior.cv=0.30, initial.yrs=5, current.yrs=5, ...) {
  eq <- eqsim(object)
  bmsy <- as.numeric(eq["bmsy"])[1]
  b0   <- as.numeric(eq["b0"])[1]

  fmsy.name <- intersect(c("fmsyMedianC", "fmsy"), dimnames(eq)[[1]])
  if (length(fmsy.name) == 0L)
    stop("eqsim(object) must contain 'fmsyMedianC' or 'fmsy'")
  fmsy <- as.numeric(eq[fmsy.name[1]])[1]

  pt <- .jabbaPellatParams(fmsy=fmsy, bmsy=bmsy, b0=b0)
  dep <- .stockDepletion(object, b0=b0, initial.yrs=initial.yrs, current.yrs=current.yrs)

  pars <- stats::setNames(
    c(
      as.numeric(pt["r"])[1],
      as.numeric(dep["psi"])[1],
      as.numeric(pt["shape"])[1],
      as.numeric(pt["k"])[1],
      as.numeric(dep["current"])[1]
    ),
    c("r", "psi", "shape", "k", "current")
  )

  .jabbaPriorObject(pars=pars, prior.cv=prior.cv)
}

.jabbaPriorsFishLife <- function(object, prior.cv=0.30, initial.yrs=5, current.yrs=5, b0.mult=8, shape=0.40, ...) {
  fl <- fishlife(object)
  r <- as.numeric(fl["r"])
  k <- b0.mult * max(as.numeric(FLCore::catch(object)), na.rm=TRUE)

  dep <- .stockDepletion(object, b0=k, initial.yrs=initial.yrs, current.yrs=current.yrs)
  pars <- c(
    .named_num(r, "r"),
    .named_num(dep["psi"], "psi"),
    .named_num(shape, "shape"),
    .named_num(k, "k"),
    .named_num(dep["current"], "current")
  )
  .jabbaPriorObject(pars=pars, prior.cv=prior.cv)
}

.jabbaPriorsFLBRP <- function(object, sr="bevholtSV", prior.cv=0.30, initial.yrs=5, current.yrs=5, ...) {
  eq <- FLRebuild::eql(object, model=sr, ...)
  rp <- FLBRP::brp(eq)

  refnames <- dimnames(FLCore::params(rp))[[1]]
  fmsy.name <- intersect(c("fmsy", "harvest", "Fmsy"), refnames)
  bmsy.name <- intersect(c("bmsy", "ssbmsy", "Bmsy"), refnames)
  b0.name   <- intersect(c("b0", "virgin", "B0"), refnames)

  if (length(fmsy.name) > 0L && length(bmsy.name) > 0L && length(b0.name) > 0L) {
    fmsy <- as.numeric(FLCore::params(rp)[fmsy.name[1]])[1]
    bmsy <- as.numeric(FLCore::params(rp)[bmsy.name[1]])[1]
    b0   <- as.numeric(FLCore::params(rp)[b0.name[1]])[1]
  } else {
    # Many FLBRP outputs store these in refpts (msy/virgin x harvest/ssb) and
    # keep params as SR coefficients (e.g. a, b).
    refs <- FLCore::refpts(rp)
    fmsy <- as.numeric(refs["msy", "harvest"])[1]
    bmsy <- as.numeric(refs["msy", "ssb"])[1]
    b0   <- as.numeric(refs["virgin", "ssb"])[1]

    if (any(is.na(c(fmsy, bmsy, b0))))
      stop("Unable to extract fmsy, bmsy, and b0-like quantities from FLBRP result")
  }

  pt <- .jabbaPellatParams(fmsy=fmsy, bmsy=bmsy, b0=b0)
  dep <- .stockDepletion(object, b0=b0, initial.yrs=initial.yrs, current.yrs=current.yrs)

  pars <- c(
    .named_num(pt["r"], "r"),
    .named_num(dep["psi"], "psi"),
    .named_num(pt["shape"], "shape"),
    .named_num(pt["k"], "k"),
    .named_num(dep["current"], "current")
  )
  .jabbaPriorObject(pars=pars, prior.cv=prior.cv)
}

#' @rdname jabbaInput
#' @exportMethod jabbaInput
setMethod("jabbaInput", signature(object="FLStock"),
function(object, index="ebiomass", catch.fun=FLCore::catch, ...) {
  catch <- data.frame(
    year  = as.integer(dimnames(FLCore::catch(object))$year),
    catch = as.numeric(c(catch.fun(object)))
  )

  if (is.character(index) && length(index) == 1L) {
    if (!exists(index, mode="function"))
      stop("Index function '", index, "' not found")
    idx.fun <- get(index, mode="function")
    idx <- .asIndexDf(idx.fun(object), value.name="index")
  } else if (is.function(index)) {
    idx <- .asIndexDf(index(object), value.name="index")
  } else {
    idx <- .asIndexDf(index, value.name="index")
  }

  list(catch=catch, index=idx)
})

.runJABBA_core <- function(catch, priors, index=NULL, model="Pella_m", assessment="", scenario="", q_bounds=NULL, sigma.est=TRUE, fixed.obsE=0.1, sigma.proc=TRUE, fixed.procE=0.3, igamma=c(3, 0.1), currentDepletion="", initialDepletion=NA, quick=TRUE, nc=3, ...) {
  .validateJabbaPriors(priors)
  .validateJabbaInput(catch=catch, index=index)

  pr    <- priors$pr
  pr.sd <- priors$pr.sd
  r <- unlist(c(pr[c("r")]))
  r.prior <- c(r, pr.sd["r"])
  psi <- unlist(c(pr[c("psi")]))
  if (is.na(psi))
    psi <- 0.9
  psi.prior <- c(psi, pr.sd["psi"])
  shape <- unlist(c(pr[c("shape")]))
  shape.cv <- pr.sd["shape"]

  k.prior <- NA
  if ("k" %in% dimnames(pr)$params) {
    k <- unlist(c(pr["k"]))
    k.prior <- c(k, pr.sd["k"])
  }

  args <- if (!is.null(q_bounds)) list(q_bounds=q_bounds) else list()
  args <- c(args, list(
    scenario   = scenario,
    assessment = assessment,
    model.type = model,
    BmsyK      = shape,
    shape.CV   = shape.cv,
    catch      = catch,
    cpue       = index,
    r.prior    = r.prior,
    K.prior    = k.prior,
    psi.prior  = psi.prior,
    sigma.proc = sigma.proc,
    sigma.est  = sigma.est,
    fixed.obsE = fixed.obsE,
    igamma     = igamma,
    verbose    = FALSE
  ))
  args <- args[!vapply(args, function(x) length(x) == 1 && all(is.na(x)), logical(1))]

  if (length(currentDepletion) > 0 && substr(currentDepletion[1], 1, 1) == "b") {
    if (!"current" %in% dimnames(pr)$params)
      stop("currentDepletion requested but prior object has no 'current'")
    args <- c(args, list(
      b.prior = c(c(pr["current"]), pr.sd["current"], max(catch$year), "bbmsy")
    ))
  }
  if (length(currentDepletion) > 0 && substr(currentDepletion[1], 1, 1) == "f") {
    if (!"current" %in% dimnames(pr)$params)
      stop("currentDepletion requested but prior object has no 'current'")
    args <- c(args, list(
      b.prior = c(c(pr["current"]), pr.sd["current"], max(catch$year), "ffmsy")
    ))
  }

  if (!is.na(initialDepletion))
    args <- c(args, list(psi.prior=c(initialDepletion, pr.sd["psi"])))

  input <- try(do.call(JABBA::build_jabba, args), silent=TRUE)
  if (inherits(input, "try-error"))
    return(NULL)

  fit <- try(JABBA::fit_jabba(input, quickmcmc=quick, verbose=FALSE, nc=nc), silent=TRUE)
  if (inherits(fit, "try-error"))
    return(list(input=input, priors=priors))

  list(input=input, fit=fit, priors=priors)
}

#' @rdname jabbaPriors
#' @exportMethod jabbaPriors
setMethod("jabbaPriors", signature(object="FLStock"),
function(object, method=c("ices", "fishlife", "flbrp"),
         prior.cv=0.30, initial.yrs=5, current.yrs=5, sr="bevholtSV",
         b0.mult=8, shape=0.40, ...) {
  method <- match.arg(method)
  switch(method,
    ices = .jabbaPriorsICES(object=object, prior.cv=prior.cv, initial.yrs=initial.yrs, current.yrs=current.yrs, ...),
    fishlife = .jabbaPriorsFishLife(object=object, prior.cv=prior.cv, initial.yrs=initial.yrs, current.yrs=current.yrs, b0.mult=b0.mult, shape=shape, ...),
    flbrp = .jabbaPriorsFLBRP(object=object, sr=sr, prior.cv=prior.cv, initial.yrs=initial.yrs, current.yrs=current.yrs, ...)
  )
})

#' @rdname runJABBA
#' @exportMethod runJABBA
setMethod("runJABBA", signature(object="FLStock"),
function(object, method=c("ices", "fishlife", "flbrp"), priors=NULL, index="ebiomass",
         model="Pella_m", assessment=as.character(FLCore::name(object)), scenario="",
         q_bounds=NULL, sigma.est=TRUE, fixed.obsE=0.1, sigma.proc=TRUE, fixed.procE=0.3,
         igamma=c(3, 0.1), currentDepletion="", initialDepletion=NA, quick=TRUE, nc=3, ...) {
  method <- match.arg(method)
  inp <- jabbaInput(object=object, index=index, ...)
  if (is.null(priors)) {
    priors <- jabbaPriors(object=object, method=method, ...)
  } else {
    .validateJabbaPriors(priors)
  }

  .runJABBA_core(
    catch=inp$catch, index=inp$index, priors=priors, model=model, assessment=assessment,
    scenario=scenario, q_bounds=q_bounds, sigma.est=sigma.est, fixed.obsE=fixed.obsE,
    sigma.proc=sigma.proc, fixed.procE=fixed.procE, igamma=igamma,
    currentDepletion=currentDepletion, initialDepletion=initialDepletion, quick=quick, nc=nc, ...
  )
})

#' @rdname runJABBA
#' @exportMethod runJABBA
setMethod("runJABBA", signature(object="data.frame"),
function(object, priors, index=NULL, model="Pella_m", assessment="", scenario="", q_bounds=NULL,
         sigma.est=TRUE, fixed.obsE=0.1, sigma.proc=TRUE, fixed.procE=0.3, igamma=c(3, 0.1),
         currentDepletion="", initialDepletion=NA, quick=TRUE, nc=3, ...) {
  .runJABBA_core(
    catch=object, index=index, priors=priors, model=model, assessment=assessment,
    scenario=scenario, q_bounds=q_bounds, sigma.est=sigma.est, fixed.obsE=fixed.obsE,
    sigma.proc=sigma.proc, fixed.procE=fixed.procE, igamma=igamma,
    currentDepletion=currentDepletion, initialDepletion=initialDepletion, quick=quick, nc=nc, ...
  )
})
