#' @rdname jabbaInput
setGeneric("jabbaInput", function(object, ...) standardGeneric("jabbaInput"))

#' JABBA prior means and standard deviations
#'
#' @param object An \code{FLStock} or \code{FLBRP}.
#' @param eq An \code{FLBRP} used with an \code{FLStock} to set \(r\), shape
#'   and \(K\) from equilibrium reference points (jabba-f).
#' @param method \code{"ices"}, \code{"fishlife"}, or \code{"flbrp"} when
#'   \code{eq} is not supplied.
#' @param stock \code{FLStock} required when \code{object} is an \code{FLBRP}.
#' @param biomass \code{"eb"} (default) or \code{"ssb"} for \(B_{\mathrm{MSY}}\)
#'   and \(B_0\).
#' @param ... Passed to the prior constructors.
#'
#' @return A list with \code{pr} and \code{pr.sd} (\code{FLPar}).
#'
#' @export
#' @rdname jabbaPriors
setGeneric("jabbaPriors",
  function(object, ...)
    standardGeneric("jabbaPriors"))

#' Fit JABBA from an FLStock and FLBRP
#'
#' @description
#' One call: catch and the biomass index come from the \code{FLStock};
#' \(F_{\mathrm{MSY}}\), \(B_{\mathrm{MSY}}\) and \(B_0\) come from the
#' \code{FLBRP}. Default index is exploitable biomass (\code{ebiomass}).
#'
#' @param object An \code{FLStock} (catch and index) or a catch \code{data.frame}.
#' @param eq An \code{FLBRP} used for production-function priors. If omitted,
#'   priors are built with \code{method}.
#' @param index Relative biomass index: a function (default \code{ebiomass}),
#'   a function name (\code{"ebiomass"}, \code{"ssb"}), an \code{FLQuant},
#'   or a \code{data.frame} with \code{year} and \code{index}.
#' @param f Optional fishing-mortality series (jabba-f auxiliary data).
#'   \code{NULL} (default) omits it. \code{ffmsy}, \code{TRUE} or
#'   \code{"ffmsy"} uses the jabba-f series
#'   \eqn{(1-e^{-F})/(1-e^{-F_{\mathrm{MSY}}})}. \code{fbar} or \code{"f"}
#'   passes absolute \(F\). An \code{FLQuant} or \code{data.frame} is used
#'   as supplied.
#' @param output What to return: \code{"jabba"} (JABBA fit; default),
#'   \code{"mpb"} (mpb \code{biodyn}), or \code{"list"} (input, fit, priors).
#' @param method Prior method when \code{eq} is not supplied. Default \code{"ices"}.
#' @param priors Optional list with \code{pr} and \code{pr.sd} (\code{FLPar}).
#' @param biomass \code{"eb"} (default) or \code{"ssb"} when deriving shape,
#'   \(K\) and depletion from \code{eq}.
#' @param model JABBA model type. Default \code{"Pella_m"}.
#' @param assessment,scenario Labels passed to JABBA.
#' @param proc.dyn If \code{TRUE}, scale process SD by depletion
#'   (\eqn{\sigma_t = \sigma \times \mathrm{relsig}(B/B_{\mathrm{MSY}})}).
#'   Requires GitHub JABBA with \code{proc.dyn}.
#' @param pe Named numeric vector of \code{\link{relsig}} coefficients
#'   when \code{proc.dyn = TRUE}, e.g. \code{pe = c(a = 0.2, c = 0.35, d = 2)}.
#'   Unnamed length 3 is taken as \code{(a, c, d)}. Defaults match JABBA
#'   (\code{a = 0.128}, \code{c = 0.223}, \code{d = 1.326}).
#' @param sigma.proc Estimate process SD (\code{TRUE}, default) or a numeric
#'   value to fix it (e.g. \code{0.07}).
#' @param igamma Inverse-gamma prior on process variance when
#'   \code{sigma.proc = TRUE}. Default \code{c(3, 0.1)}.
#' @param quick If \code{TRUE}, use JABBA's short MCMC.
#' @param nc Number of MCMC chains.
#' @param ... Passed to \code{jabbaInput()}, \code{jabbaPriors()} and JABBA.
#'
#' @return A JABBA fit (default), an mpb \code{biodyn} if \code{output = "mpb"},
#'   or a list with \code{input}, \code{fit} and \code{priors}.
#'
#' @examples
#' \dontrun{
#' data(ple4, package = "FLCore")
#' data(ple4brp, package = "FLBRP")
#' fit  <- runJABBA(ple4, ple4brp, index = ebiomass)
#' fitF <- runJABBA(ple4, ple4brp, index = ebiomass, f = ffmsy)
#' bd   <- runJABBA(ple4, ple4brp, index = ebiomass, output = "mpb")
#'
#' # Process error: estimate SD (default), fix it, or scale by B/BMSY
#' fitFix <- runJABBA(ple4, ple4brp, index = ebiomass, sigma.proc = 0.07)
#' fitPE  <- runJABBA(ple4, ple4brp, index = ebiomass, proc.dyn = TRUE)
#' fitPE2 <- runJABBA(ple4, ple4brp, index = ebiomass, proc.dyn = TRUE,
#'                    pe = c(a = 0.2, c = 0.35, d = 2))
#' bdPE   <- jabba2biodyn(fitPE)
#' plotPe(jabbaTs(fitPE))
#' }
#'
#' @seealso \code{\link{jabbaInput}}, \code{\link{jabbaPriors}},
#'   \code{\link{compareJabbaPE}}, \code{\link{jabba2biodyn}}
#' @export
#' @rdname runJABBA
setGeneric("runJABBA", function(object, ...) standardGeneric("runJABBA"))

.validateJabbaPriors<-function(priors) {
  if (!is.list(priors) || !all(c("pr", "pr.sd") %in% names(priors)))
    stop("priors must be a list with elements 'pr' and 'pr.sd'")

  if (!methods::is(priors$pr, "FLPar"))
    stop("priors$pr must be an FLPar")

  if (!methods::is(priors$pr.sd, "FLPar"))
    stop("priors$pr.sd must be an FLPar")

  req=c("r", "psi", "shape")
  pnm=dimnames(priors$pr)$params
  snm=dimnames(priors$pr.sd)$params

  if (!all(req %in% pnm))
    stop("priors$pr must contain parameters: ", paste(req, collapse=", "))
  if (!all(req %in% snm))
    stop("priors$pr.sd must contain parameters: ", paste(req, collapse=", "))

  opt=intersect(c("k", "current"), pnm)
  if (length(opt) > 0L && !all(opt %in% snm))
    stop("priors$pr.sd must also contain optional parameters found in priors$pr: ",
         paste(opt, collapse=", "))

  invisible(priors)
}

.validateJabbaInput<-function(catch, index=NULL) {
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

.stockDepletion<-function(object, b0, initial.yrs=5, current.yrs=5, biomass="e") {
  ts=tseries(object)
  yrs0=seq(min(ts$year), length.out=initial.yrs)
  yrs1=seq(max(ts$year) - current.yrs + 1, max(ts$year))
  
  if (substr(biomass,1,1)=="e")
  c(psi     = mean(ts$eb[ts$year %in% yrs0], na.rm=TRUE) / b0,
    current = mean(ts$eb[ts$year %in% yrs1], na.rm=TRUE) / b0)
  else
    c(psi     = mean(ts$ssb[ts$year %in% yrs0], na.rm=TRUE) / b0,
      current = mean(ts$ssb[ts$year %in% yrs1], na.rm=TRUE) / b0)
}

.index2df<-function(x, value.name="index") {
  if (methods::is(x, "FLQuant")) {
    yrs =as.integer(dimnames(x)$year)
    vals=as.numeric(c(x))
    return(data.frame(year=yrs, index=vals))
  }
  if (is.data.frame(x)) {
    if (!all(c("year", value.name) %in% names(x)))
      stop("index data.frame must contain columns 'year' and '", value.name, "'")
    return(x[, c("year", value.name), drop=FALSE])
  }
  stop("Unsupported index object passed to .index2df()")
}

.asAuxDf<-function(x, value.name="f") {
  if (methods::is(x, "FLQuant")) {
    yrs =as.integer(dimnames(x)$year)
    vals=as.numeric(c(x))
    out=data.frame(year=yrs)
    out[[value.name]]=vals
    return(out)
  }
  if (is.data.frame(x)) {
    if ("year" %in% names(x) && value.name %in% names(x))
      return(x[, c("year", value.name), drop=FALSE])
    if ("year" %in% names(x) && ncol(x) >= 2L) {
      out=x[, 1:2, drop=FALSE]
      names(out)=c("year", value.name)
      return(out)
    }
    stop("F data.frame must contain 'year' and '", value.name, "'")
  }
  stop("Unsupported F object")
}

.jabbaFmsy<-function(object, eq=NULL) {
  if (!is.null(eq) && methods::is(eq, "FLBRP")) {
    v=as.numeric(FLBRP::refpts(eq)["msy", "harvest"])[1]
    if (is.finite(v) && v > 0)
      return(v)
  }
  bm=try(benchmark(object), silent=TRUE)
  if (!inherits(bm, "try-error") && !is.null(bm) &&
      "fmsy" %in% dimnames(bm)$params) {
    v=as.numeric(bm["fmsy"])[1]
    if (is.finite(v) && v > 0)
      return(v)
  }
  stop("ffmsy needs FMSY from eq or benchmark(object)")
}

#' Harvest rate relative to \(F_{\mathrm{MSY}}\) (jabba-f)
#'
#' \deqn{(1-e^{-F})/(1-e^{-F_{\mathrm{MSY}}})}
#' with \(F\) from \code{fbar}. \(F_{\mathrm{MSY}}\) comes from \code{eq}
#' when supplied, otherwise from \code{benchmark(object)}.
#'
#' @param object An \code{FLStock}.
#' @param eq Optional \code{FLBRP}.
#' @param fmsy Optional scalar; overrides \code{eq} and \code{benchmark}.
#' @param ... Unused.
#'
#' @return An \code{FLQuant}.
#'
#' @examples
#' \dontrun{
#' ffmsy(ple4, eq = ple4brp)
#' runJABBA(ple4, ple4brp, index = ebiomass, f = ffmsy)
#' }
#'
#' @seealso \code{\link{runJABBA}}, \code{\link{ebiomass}}, \code{\link{fbar}}
#' @export
#' @rdname ffmsy
setGeneric("ffmsy", function(object, ...) standardGeneric("ffmsy"))

#' @rdname ffmsy
#' @export
setMethod("ffmsy", signature(object="FLStock"),
function(object, eq=NULL, fmsy=NULL, ...) {
  if (is.null(fmsy))
    fmsy=.jabbaFmsy(object, eq)
  fv=FLCore::fbar(object)
  (1 - exp(-fv)) / (1 - exp(-fmsy))
})

.jabbaAuxiliary<-function(object, f=NULL, eq=NULL) {
  if (is.null(f) || isFALSE(f))
    return(NULL)

  if (isTRUE(f) || (is.character(f) && length(f) == 1L &&
                    tolower(f) %in% c("ffmsy", "f/fmsy")) ||
      (is.function(f) && identical(f, ffmsy))) {
    aux=.asAuxDf(ffmsy(object, eq=eq), value.name="ffmsy")
    return(list(auxiliary=aux, auxiliary.type="ffmsy"))
  }

  if (is.character(f) && length(f) == 1L && tolower(f) %in% c("f", "fbar")) {
    aux=.asAuxDf(FLCore::fbar(object), value.name="f")
    return(list(auxiliary=aux, auxiliary.type="f"))
  }

  if (is.function(f)) {
    aux=.asAuxDf(f(object), value.name="f")
    type=if (identical(f, FLCore::fbar) || identical(f, FLCore::harvest))
      "f" else "ffmsy"
    if (identical(type, "ffmsy"))
      names(aux)[2]="ffmsy"
    return(list(auxiliary=aux, auxiliary.type=type))
  }

  if (is.data.frame(f)) {
    type=if ("ffmsy" %in% names(f)) "ffmsy" else "f"
    aux=.asAuxDf(f, value.name=names(f)[names(f) != "year"][1])
    names(aux)[2]=type
    return(list(auxiliary=aux, auxiliary.type=type))
  }

  aux=.asAuxDf(f, value.name="f")
  list(auxiliary=aux, auxiliary.type="f")
}

.asIndexDf<-function(x, value.name="index") {
  if (exists(".index2df", mode="function", inherits=TRUE))
    return(.index2df(x, value.name=value.name))

  if (methods::is(x, "FLQuant")) {
    yrs =as.integer(dimnames(x)$year)
    vals=as.numeric(c(x))
    return(data.frame(year=yrs, index=vals))
  }
  if (is.data.frame(x)) {
    if (!all(c("year", value.name) %in% names(x)))
      stop("index data.frame must contain columns 'year' and '", value.name, "'")
    return(x[, c("year", value.name), drop=FALSE])
  }
  stop("Unsupported index object passed to jabbaInput(index=...)")
}

.jabbaPellatParams<-function(fmsy, bmsy, b0, interval=NULL) {
  fmsy=as.numeric(fmsy)[1]
  bmsy=as.numeric(bmsy)[1]
  b0  =as.numeric(b0)[1]
  shape=bmsy / b0

  if (is.null(interval))
    interval=if (shape < 0.37) c(0.10, 0.99) else c(1.01, 10.0)

  m=stats::optimize(
    function(x) abs(shape - (1 / x)^(1 / (x - 1))),
    interval=interval
  )$minimum
  r=fmsy * (m - 1) / (1 - 1 / m)
  FLCore::FLPar(c(r=r, shape=shape, k=b0))
}

.jabbaPriorObject<-function(pars, sds=NULL, prior.cv=0.30) {
  if (is.list(pars))
    pars=unlist(pars, use.names=TRUE)
  pars=stats::setNames(as.numeric(pars), names(pars))

  req=c("r", "psi", "shape")
  if (!all(req %in% names(pars)))
    stop(
      "Prior means must contain: ", paste(req, collapse=", "),
      ". Found: ", paste(names(pars), collapse=", ")
    )

  keep=intersect(c("r", "psi", "shape", "k", "current"), names(pars))
  pr=FLCore::FLPar(pars[keep])

  if (is.null(sds)) {
    pr.sd=pr %/% pr * prior.cv
  } else {
    if (is.list(sds))
      sds=unlist(sds, use.names=TRUE)
    sds=stats::setNames(as.numeric(sds), names(sds))
    if (!all(keep %in% names(sds)))
      stop(
        "Prior sds must contain: ", paste(keep, collapse=", "),
        ". Found: ", paste(names(sds), collapse=", ")
      )
    pr.sd=FLCore::FLPar(sds[keep])
  }

  list(pr=pr, pr.sd=pr.sd)
}

.named_num<-function(x, nm) {
  stats::setNames(as.numeric(x)[1], nm)
}

.jabbaPriorsICES<-function(object, prior.cv=0.30, initial.yrs=5, current.yrs=5, ...) {
  eq=eqsim(object)
  bmsy=as.numeric(eq["bmsy"])[1]
  b0  =as.numeric(eq["b0"])[1]

  fmsy.name=intersect(c("fmsyMedianC", "fmsy"), dimnames(eq)[[1]])
  if (length(fmsy.name) == 0L)
    stop("eqsim(object) must contain 'fmsyMedianC' or 'fmsy'")
  fmsy=as.numeric(eq[fmsy.name[1]])[1]

  pt=.jabbaPellatParams(fmsy=fmsy, bmsy=bmsy, b0=b0)
  dep=.stockDepletion(object, b0=b0, initial.yrs=initial.yrs, current.yrs=current.yrs)

  pars=stats::setNames(
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

.jabbaPriorsFishLife<-function(object, prior.cv=0.30, initial.yrs=5, current.yrs=5, b0.mult=8, shape=0.40, ...) {
  fl=fishlife(object)
  r=as.numeric(fl["r"])
  k=b0.mult * max(as.numeric(FLCore::catch(object)), na.rm=TRUE)

  dep=.stockDepletion(object, b0=k, initial.yrs=initial.yrs, current.yrs=current.yrs)
  pars=c(
    .named_num(r, "r"),
    .named_num(dep["psi"], "psi"),
    .named_num(shape, "shape"),
    .named_num(k, "k"),
    .named_num(dep["current"], "current")
  )
  .jabbaPriorObject(pars=pars, prior.cv=prior.cv)
}

.jabbaPriorsFromBRP<-function(eq, stock, prior.cv=0.30, initial.yrs=5,
                             current.yrs=5, biomass="eb", ...) {
  if (!methods::is(eq, "FLBRP"))
    stop("'eq' must be an FLBRP")
  if (!methods::is(stock, "FLStock"))
    stop("'stock' must be an FLStock")

  rp=try(FLCore::refpts(eq)[c("msy", "virgin"), c("harvest", "ssb")], silent=TRUE)
  if (inherits(rp, "try-error") || any(!is.finite(as.numeric(rp))))
    eq=FLBRP::brp(eq)

  if (substr(biomass, 1, 1) == "e") {
    refs=refptsEB(eq)
    bmsy=as.numeric(refs["msy",    "eb"])[1]
    b0  =as.numeric(refs["virgin", "eb"])[1]
  } else {
    refs=FLCore::refpts(eq)
    bmsy=as.numeric(refs["msy",    "ssb"])[1]
    b0  =as.numeric(refs["virgin", "ssb"])[1]
  }
  fmsy=as.numeric(refs["msy", "harvest"])[1]
  if (any(!is.finite(c(fmsy, bmsy, b0))) || any(c(fmsy, bmsy, b0) <= 0))
    stop("Unable to extract finite positive fmsy, bmsy and b0 from FLBRP")

  pt=.jabbaPellatParams(fmsy=fmsy, bmsy=bmsy, b0=b0)
  dep=.stockDepletion(stock, b0=b0, initial.yrs=initial.yrs,
                      current.yrs=current.yrs, biomass=biomass)

  pars=c(
    .named_num(pt["r"], "r"),
    .named_num(dep["psi"], "psi"),
    .named_num(pt["shape"], "shape"),
    .named_num(pt["k"], "k"),
    .named_num(dep["current"], "current")
  )
  .jabbaPriorObject(pars=pars, prior.cv=prior.cv)
}

.jabbaPriorsFLBRP<-function(object, sr="bevholtSV", prior.cv=0.30, initial.yrs=5,
                            current.yrs=5, biomass="eb", ...) {
  eq=FLRebuild::eql(object, model=sr, ...)
  .jabbaPriorsFromBRP(eq=eq, stock=object, prior.cv=prior.cv,
                      initial.yrs=initial.yrs, current.yrs=current.yrs,
                      biomass=biomass, ...)
}

#' @rdname jabbaInput
#' @exportMethod jabbaInput
setMethod("jabbaInput", signature(object="FLStock"),
function(object, index=ebiomass, f=NULL, eq=NULL, catch.fun=FLCore::catch, ...) {
  catch=data.frame(
    year  = as.integer(dimnames(FLCore::catch(object))$year),
    catch = as.numeric(c(catch.fun(object)))
  )

  if (is.character(index) && length(index) == 1L) {
    if (!exists(index, mode="function"))
      stop("Index function '", index, "' not found")
    idx.fun=get(index, mode="function")
    idx=.asIndexDf(idx.fun(object), value.name="index")
  } else if (is.function(index)) {
    idx=.asIndexDf(index(object), value.name="index")
  } else {
    idx=.asIndexDf(index, value.name="index")
  }

  aux=.jabbaAuxiliary(object, f=f, eq=eq)
  out=list(catch=catch, index=idx)
  if (!is.null(aux)) {
    out$auxiliary=aux$auxiliary
    out$auxiliary.type=aux$auxiliary.type
  }
  out
})

.jabbaBuildFormals<-function() {
  if (!requireNamespace("JABBA", quietly=TRUE))
    stop("Package 'JABBA' is required. Install with remotes::install_github(\"jabbamodel/JABBA\")")
  names(formals(JABBA::build_jabba))
}

.jabbaHasProcDyn<-function() {
  "proc.dyn" %in% .jabbaBuildFormals()
}

.jabbaRelsigDefault <- c(a = 0.128, c = 0.223, d = 1.326)

.jabbaParsePeVector<-function(pe) {
  if (is.null(pe))
    return(NULL)
  if (is.list(pe))
    pe=unlist(pe)
  if (!is.numeric(pe) || !length(pe))
    stop("'pe' must be numeric, e.g. pe=c(a=0.2, c=0.35, d=2)")
  nms=names(pe)
  out=as.list(.jabbaRelsigDefault)
  if (is.null(nms) || !any(nzchar(nms))) {
    if (length(pe) != 3L)
      stop("'pe' must be named (a, c, d) or a length-3 vector in that order")
    out$a=unname(pe[1])
    out$c=unname(pe[2])
    out$d=unname(pe[3])
    return(out)
  }
  nms=tolower(trimws(nms))
  bad=setdiff(nms[nzchar(nms)], c("a", "c", "d"))
  if (length(bad))
    stop("'pe' names must be a, c, d, e.g. pe=c(a=0.2, c=0.35, d=2)")
  names(pe)=nms
  for (nm in c("a", "c", "d")) {
    if (nm %in% nms)
      out[[nm]]=unname(pe[nms == nm][1])
  }
  out
}

.jabbaRelsigPars<-function(pe=NULL, ...) {
  dots=list(...)
  if (is.null(pe) && !is.null(dots$pe))
    pe=dots$pe
  vec=.jabbaParsePeVector(pe)
  if (is.null(vec))
    vec=as.list(.jabbaRelsigDefault)
  a=as.numeric(vec$a)[1]
  cc=as.numeric(vec$c)[1]
  d=as.numeric(vec$d)[1]
  if (!is.finite(a) || a <= 0 || a >= 1)
    stop("pe['a'] (relsig floor) must be in (0, 1)")
  if (!is.finite(cc) || cc <= 0)
    stop("pe['c'] (relsig inflection) must be > 0")
  if (!is.finite(d) || d <= 0)
    stop("pe['d'] (relsig steepness) must be > 0")
  list(a=a, c=cc, d=d)
}

.jabbaDotsWithoutRelsig<-function(dots) {
  nms=names(dots)
  if (is.null(nms))
    return(dots)
  dots[nms != "pe"]
}

.jabbaRelsigIsDefault<-function(rs) {
  isTRUE(all.equal(c(rs$a, rs$c, rs$d), unname(.jabbaRelsigDefault),
                   tolerance=1e-8))
}

.patchJabbaJagsRelsig<-function(file, a, c, d) {
  if (!file.exists(file))
    return(FALSE)
  txt=readLines(file, warn=FALSE)
  if (!any(grepl("a_pe", txt)))
    return(FALSE)
  txt=sub("a_pe\\s*<-\\s*[0-9.]+", sprintf("a_pe <- %s", a), txt)
  txt=sub("c_pe\\s*<-\\s*[0-9.]+", sprintf("c_pe <- %s", c), txt)
  txt=sub("d_pe\\s*<-\\s*[0-9.]+", sprintf("d_pe <- %s", d), txt)
  writeLines(txt, file)
  TRUE
}

.fitJabba<-function(input, quick, nc, silent, rs, proc.dyn) {
  fitfun=function() {
    try(JABBA::fit_jabba(input, quickmcmc=quick, verbose=FALSE, nc=nc),
        silent=silent)
  }
  if (!isTRUE(proc.dyn) || .jabbaRelsigIsDefault(rs))
    return(fitfun())
  ns=asNamespace("JABBA")
  if (!exists("jabba2jags", envir=ns, inherits=FALSE)) {
    warning("Cannot set pe: JABBA::jabba2jags not found")
    return(fitfun())
  }
  orig=get("jabba2jags", envir=ns, inherits=FALSE)
  patched=function(jbinput, dir) {
    orig(jbinput, dir)
    ok=.patchJabbaJagsRelsig(file.path(dir, "JABBA.jags"), rs$a, rs$c, rs$d)
    if (!ok)
      warning("Could not write pe into JABBA.jags")
    invisible(NULL)
  }
  unlockBinding("jabba2jags", ns)
  assign("jabba2jags", patched, envir=ns)
  lockBinding("jabba2jags", ns)
  on.exit({
    unlockBinding("jabba2jags", ns)
    assign("jabba2jags", orig, envir=ns)
    lockBinding("jabba2jags", ns)
  }, add=TRUE)
  fitfun()
}

.runJABBA_core<-function(catch, priors, index=NULL, model="Pella_m", assessment="", scenario="", q_bounds=NULL, sigma.est=TRUE, fixed.obsE=0.1, sigma.proc=TRUE, fixed.procE=0.3, igamma=c(3, 0.1), currentDepletion="", initialDepletion=NA, proc.dyn=FALSE, pe=NULL, auxiliary=NULL, auxiliary.type="ffmsy", quick=TRUE, nc=3, silent=TRUE,...) {
  .validateJabbaPriors(priors)
  .validateJabbaInput(catch=catch, index=index)

  pr   =priors$pr
  pr.sd=priors$pr.sd
  r=unlist(c(pr[c("r")]))
  r.prior=c(r, pr.sd["r"])
  psi=unlist(c(pr[c("psi")]))
  if (is.na(psi))
    psi=0.9
  psi.prior=c(psi, pr.sd["psi"])
  shape=unlist(c(pr[c("shape")]))
  shape.cv=pr.sd["shape"]

  k.prior=NA
  if ("k" %in% dimnames(pr)$params) {
    k=unlist(c(pr["k"]))
    k.prior=c(k, pr.sd["k"])
  }

  args=if (!is.null(q_bounds)) list(q_bounds=q_bounds) else list()
  args=c(args, list(
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
  args=args[!vapply(args, function(x) length(x) == 1 && all(is.na(x)), logical(1))]

  jb.formals=.jabbaBuildFormals()
  if (isTRUE(proc.dyn) && !"proc.dyn" %in% jb.formals)
    stop("Installed JABBA has no proc.dyn argument. Update with ",
         "remotes::install_github(\"jabbamodel/JABBA\")")
  if ("proc.dyn" %in% jb.formals)
    args$proc.dyn=isTRUE(proc.dyn)

  if (!is.null(auxiliary)) {
    args$auxiliary=auxiliary
    args$auxiliary.type=auxiliary.type
    args$auxiliary.sigma=TRUE
    args$auxiliary.obsE=0.3
    args$auxiliary.lag=0
  }

  extra=list(...)
  rs=.jabbaRelsigPars(pe=pe)
  extra=.jabbaDotsWithoutRelsig(extra)
  extra=extra[setdiff(names(extra), names(args))]
  extra=extra[names(extra) %in% jb.formals]
  if (length(extra))
    args=c(args, extra)

  if (!isTRUE(proc.dyn) && !.jabbaRelsigIsDefault(rs))
    warning("'pe' is ignored unless proc.dyn = TRUE")

  for (nm in c("a_pe", "c_pe", "d_pe")) {
    if (nm %in% jb.formals) {
      key=switch(nm, a_pe="a", c_pe="c", d_pe="d")
      args[[nm]]=rs[[key]]
    }
  }

  if (length(currentDepletion) > 0 && substr(currentDepletion[1], 1, 1) == "b") {
    if (!"current" %in% dimnames(pr)$params)
      stop("currentDepletion requested but prior object has no 'current'")
    args=c(args, list(
      b.prior = c(c(pr["current"]), pr.sd["current"], max(catch$year), "bbmsy")
    ))
  }
  if (length(currentDepletion) > 0 && substr(currentDepletion[1], 1, 1) == "f") {
    if (!"current" %in% dimnames(pr)$params)
      stop("currentDepletion requested but prior object has no 'current'")
    args=c(args, list(
      b.prior = c(c(pr["current"]), pr.sd["current"], max(catch$year), "ffmsy")
    ))
  }

  if (!is.na(initialDepletion))
    args=c(args, list(psi.prior=c(initialDepletion, pr.sd["psi"])))

  input=try(do.call(JABBA::build_jabba, args), silent=silent)
  if (inherits(input, "try-error"))
    return(NULL)
  if (is.list(input) && is.list(input$settings))
    input$settings$relsig=c(a=rs$a, c=rs$c, d=rs$d)

  fit=.fitJabba(input, quick=quick, nc=nc, silent=silent, rs=rs, proc.dyn=proc.dyn)
  if (inherits(fit, "try-error"))
    return(list(input=input, priors=priors))
  if (is.list(fit))
    fit$relsig=c(a=rs$a, c=rs$c, d=rs$d)

  list(input=input, fit=fit, priors=priors)
}

#' Coerce a JABBA fit to an mpb \code{biodyn}
#'
#' Copies catch, biomass, index, median \code{r}, \(K\), shape and
#' initial depletion from a JABBA fit into an mpb \code{biodyn}.
#' Accepts a JABBA fit or a wrapper with \code{fit$timeseries}.
#'
#' @param object A JABBA fit (or list with a \code{fit} element).
#'
#' @return An mpb \code{biodyn}.
#'
#' @examples
#' \dontrun{
#' fit <- runJABBA(ple4, ple4brp, index = ebiomass)
#' bd  <- jabba2biodyn(fit)
#' }
#'
#' @seealso \code{\link{runJABBA}}
#' @export
jabba2biodyn<-function(object) {
  if (!requireNamespace("mpb", quietly=TRUE))
    stop("Package 'mpb' is required. Install mpb and retry")
  if (is.null(object))
    stop("'object' cannot be NULL")
  if (is.null(object$timeseries) && !is.null(object$fit$timeseries))
    object <- object$fit
  if (is.null(object$timeseries))
    stop("Cannot coerce to biodyn: JABBA fit has no timeseries")

  res=mpb::biodyn()
  res@name=paste(object$assessment, object$scenario)
  res@desc="coerced from JABBA"

  catch.df=object$inputseries$catch
  res@catch=FLCore::as.FLQuant(data.frame(
    year=catch.df$year, data=catch.df$catch))
  res@stock=FLCore::as.FLQuant(data.frame(
    year=as.numeric(dimnames(object$timeseries)[[1]]),
    data=as.numeric(object$timeseries[, "mu", "B"])))

  if (!is.null(object$inputseries$cpue)) {
    cpue=object$inputseries$cpue
    idx=list()
    for (nm in names(cpue)[-1])
      idx[[nm]]=FLCore::as.FLQuant(data.frame(year=cpue[[1]], data=cpue[[nm]]))
    if (length(idx))
      res@indices=do.call(FLCore::FLQuants, idx)
  }

  if (!is.null(object$pars)) {
    p=object$pars
    col=if ("Median" %in% colnames(p)) "Median" else 1
    getp=function(nm) {
      if (is.null(rownames(p)) || !nm %in% rownames(p))
        return(NA_real_)
      as.numeric(p[nm, col])
    }
    res@params["r"] =getp("r")
    res@params["k"] =getp("K")
    res@params["p"] =getp("m") - 1
    res@params["b0"]=getp("psi")
  }

  if (!is.null(object$kobe))
    res@kobe=object$kobe
  if (!is.null(object$pars_posterior))
    res@posteriors=object$pars_posterior

  res
}

.jabbaFormatOutput<-function(out, output=c("jabba", "mpb", "list")) {
  output=match.arg(output)
  if (is.null(out))
    return(NULL)
  if (identical(output, "list"))
    return(out)
  fit=out$fit
  if (is.null(fit)) {
    warning("JABBA MCMC did not return a fit")
    return(out)
  }
  fit$priors=out$priors
  if (identical(output, "jabba"))
    return(fit)
  jabba2biodyn(fit)
}

#' @rdname jabbaPriors
#' @exportMethod jabbaPriors
setMethod("jabbaPriors", signature(object="FLStock"),
function(object, eq=NULL, method=c("ices", "fishlife", "flbrp"),
         prior.cv=0.30, initial.yrs=5, current.yrs=5, sr="bevholtSV",
         b0.mult=8, shape=0.40, biomass="eb", ...) {
  if (!is.null(eq)) {
    if (!methods::is(eq, "FLBRP"))
      stop("'eq' must be an FLBRP")
    return(.jabbaPriorsFromBRP(eq=eq, stock=object, prior.cv=prior.cv,
                               initial.yrs=initial.yrs, current.yrs=current.yrs,
                               biomass=biomass, ...))
  }
  method=match.arg(method)
  switch(method,
    ices = .jabbaPriorsICES(object=object, prior.cv=prior.cv, initial.yrs=initial.yrs, current.yrs=current.yrs, ...),
    fishlife = .jabbaPriorsFishLife(object=object, prior.cv=prior.cv, initial.yrs=initial.yrs, current.yrs=current.yrs, b0.mult=b0.mult, shape=shape, ...),
    flbrp = .jabbaPriorsFLBRP(object=object, sr=sr, prior.cv=prior.cv, initial.yrs=initial.yrs, current.yrs=current.yrs, biomass=biomass, ...)
  )
})

#' @rdname jabbaPriors
#' @exportMethod jabbaPriors
setMethod("jabbaPriors", signature(object="FLBRP"),
function(object, stock, prior.cv=0.30, initial.yrs=5, current.yrs=5,
         biomass="eb", ...) {
  if (missing(stock))
    stop("jabbaPriors(FLBRP) needs a stock= FLStock for depletion (psi, current)")
  .jabbaPriorsFromBRP(eq=object, stock=stock, prior.cv=prior.cv,
                      initial.yrs=initial.yrs, current.yrs=current.yrs,
                      biomass=biomass, ...)
})

#' @rdname runJABBA
#' @exportMethod runJABBA
setMethod("runJABBA", signature(object="FLStock"),
function(object, eq=NULL, index=ebiomass, f=NULL, output=c("jabba", "mpb", "list"),
         method=c("ices", "fishlife", "flbrp"), priors=NULL, biomass="eb",
         model="Pella_m", assessment=as.character(FLCore::name(object)), scenario="",
         q_bounds=NULL, sigma.est=TRUE, fixed.obsE=0.1, sigma.proc=TRUE, fixed.procE=0.3,
         igamma=c(3, 0.1), currentDepletion="", initialDepletion=NA, proc.dyn=FALSE,
         pe=NULL, quick=TRUE, nc=3, ...) {
  output=match.arg(output)
  if (!is.null(eq) && !methods::is(eq, "FLBRP"))
    stop("'eq' must be an FLBRP")
  inp=jabbaInput(object=object, index=index, f=f, eq=eq, ...)
  if (is.null(priors)) {
    if (!is.null(eq)) {
      priors=jabbaPriors(object, eq=eq, biomass=biomass, ...)
    } else {
      method=match.arg(method)
      priors=jabbaPriors(object=object, method=method, biomass=biomass, ...)
    }
  } else {
    .validateJabbaPriors(priors)
  }

  out=.runJABBA_core(
    catch=inp$catch, index=inp$index, priors=priors, model=model, assessment=assessment,
    scenario=scenario, q_bounds=q_bounds, sigma.est=sigma.est, fixed.obsE=fixed.obsE,
    sigma.proc=sigma.proc, fixed.procE=fixed.procE, igamma=igamma,
    currentDepletion=currentDepletion, initialDepletion=initialDepletion,
    proc.dyn=proc.dyn, pe=pe,
    auxiliary=inp$auxiliary, auxiliary.type=if (is.null(inp$auxiliary.type)) "ffmsy" else inp$auxiliary.type,
    quick=quick, nc=nc, ...
  )
  .jabbaFormatOutput(out, output)
})

#' @rdname runJABBA
#' @exportMethod runJABBA
setMethod("runJABBA", signature(object="data.frame"),
function(object, priors, index=NULL, f=NULL, output=c("jabba", "mpb", "list"),
         model="Pella_m", assessment="", scenario="", q_bounds=NULL,
         sigma.est=TRUE, fixed.obsE=0.1, sigma.proc=TRUE, fixed.procE=0.3, igamma=c(3, 0.1),
         currentDepletion="", initialDepletion=NA, proc.dyn=FALSE,
         pe=NULL, quick=TRUE, nc=3, ...) {
  output=match.arg(output)
  aux=NULL
  atype="ffmsy"
  if (!is.null(f)) {
    if (is.list(f) && !is.null(f$auxiliary)) {
      aux=f$auxiliary
      atype=f$auxiliary.type
    } else if (is.data.frame(f)) {
      atype=if ("ffmsy" %in% names(f)) "ffmsy" else "f"
      aux=.asAuxDf(f, value.name=names(f)[names(f) != "year"][1])
      names(aux)[2]=atype
    } else {
      stop("For a catch data.frame, 'f' must be a data.frame of F or F/FMSY")
    }
  }
  out=.runJABBA_core(
    catch=object, index=index, priors=priors, model=model, assessment=assessment,
    scenario=scenario, q_bounds=q_bounds, sigma.est=sigma.est, fixed.obsE=fixed.obsE,
    sigma.proc=sigma.proc, fixed.procE=fixed.procE, igamma=igamma,
    currentDepletion=currentDepletion, initialDepletion=initialDepletion,
    proc.dyn=proc.dyn, pe=pe, auxiliary=aux, auxiliary.type=atype, quick=quick, nc=nc, ...
  )
  .jabbaFormatOutput(out, output)
})
