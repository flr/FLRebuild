#' Fit FLICC Models from FLR Stocks
#'
#' Generic interface to prepare FLICC inputs from FLR stock objects and fit
#' \code{FLicc::fiticc()}.
#'
#' @param object Stock object.
#' @param ... Additional method arguments.
#'
#' @return A list with estimated parameters, model input and fit output.
#' @export
setGeneric("fitFLICC", function(object, ...) standardGeneric("fitFLICC"))

.vbFn <- function(age, linf, k, t0) {
  linf * (1 - exp(-k * (age - t0)))
}

.vonbFromStock <- function(stk) {
  ln <- FLCore::wt2len(FLCore::stock.wt(stk), FLCore::FLPar(a=0.1, b=3))
  ln <- as.data.frame(ln, drop=TRUE)
  ln <- ln[is.finite(ln$data), ]
  ln$age <- as.numeric(as.character(ln$age))

  fit <- stats::nls(
    data ~ .vbFn(age, linf, k, t0),
    data=ln,
    start=list(linf=max(ln$data, na.rm=TRUE), k=0.2, t0=-0.2)
  )
  stats::coef(fit)
}

.mat50FromStock <- function(stk) {
  md <- as.data.frame(FLCore::mat(stk))
  md <- md[is.finite(md$data), ]
  md$age <- as.numeric(as.character(md$age))

  byYear <- split(md, md$year)
  vapply(byYear, function(df) {
    df <- df[order(df$age), ]
    if (any(abs(df$data - 0.5) < 1e-8))
      return(df$age[which.min(abs(df$data - 0.5))])
    i <- which(df$data > 0.5)[1]
    if (is.na(i) || i == 1)
      return(NA_real_)
    x1 <- df$age[i - 1]; y1 <- df$data[i - 1]
    x2 <- df$age[i];     y2 <- df$data[i]
    x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)
  }, numeric(1))
}

.lorFromStock <- function(stk) {
  mdat <- as.data.frame(FLCore::m(stk))
  wdat <- as.data.frame(FLCore::stock.wt(stk))
  df <- merge(
    mdat, wdat,
    by=c("year", "age", "unit", "season", "area", "iter"),
    suffixes=c(".m", ".w")
  )
  df <- df[is.finite(df$data.m) & is.finite(df$data.w) & df$data.w > 0, ]
  fit <- stats::glm(log(data.m) ~ log(data.w), data=df, family=stats::gaussian())
  cf <- stats::coef(fit)
  stats::setNames(c(cf[[1]], cf[[2]]), c("m1", "m2"))
}

.lorToGis <- function(par) {
  c3 <- 1.44
  c4 <- 1
  b <- as.numeric(par["m1"])
  a <- as.numeric(par["m2"]) - c3 * log(as.numeric(par["linf"])) - c4 * log(as.numeric(par["k"]))
  FLCore::FLPar(c(m1=a, m2=b, m3=c3, m4=c4))
}

.wlFromStock <- function(stk, linf, k, t0) {
  wt <- as.numeric(FLCore::stock.wt(stk)[, 1, 1, 1, 1, 1])
  age <- as.numeric(dimnames(FLCore::stock.wt(stk))$age)
  len <- .vbFn(age, linf=linf, k=k, t0=t0)
  ok <- is.finite(wt) & is.finite(len) & wt > 0 & len > 0
  fit <- stats::lm(log(wt[ok]) ~ log(len[ok]))
  cf <- stats::coef(fit)
  c(a=exp(cf[[1]]), b=cf[[2]])
}

#' @rdname fitFLICC
#' @param lengthData Optional length-frequency FLQuant. If NULL, the method
#'   tries to build one via \code{icesdata::lenSamples()}.
#' @param nLengths Number of simulated fish for \code{lenSamples()}.
#' @param catch Total catch value used in FLICC input.
#' @param mk Initial \code{Mk} parameters passed to FLICC.
#' @param selFun Selectivity function name.
#' @param refLength Reference length passed to FLICC.
#' @param compile Logical passed to \code{FLicc::fiticc()}.
#' @param modelName Label for the FLICC model.
#' @exportMethod fitFLICC
setMethod("fitFLICC", signature(object="FLStock"),
function(object,
         lengthData=NULL,
         nLengths=100,
         catch=12000,
         mk=c(2, 0.01),
         selFun="logistic",
         refLength=5,
         compile=TRUE,
         modelName="Base: Mk ~ Length-inverse",
         ...) {
  dots <- list(...)
  if (is.null(lengthData) && !is.null(dots$length_data))
    lengthData <- dots$length_data
  if (!is.null(dots$n_lengths))
    nLengths <- dots$n_lengths
  if (!is.null(dots$sel_fun))
    selFun <- dots$sel_fun
  if (!is.null(dots$ref_length))
    refLength <- dots$ref_length
  if (!is.null(dots$model_name))
    modelName <- dots$model_name

  if (!requireNamespace("FLicc", quietly=TRUE))
    stop("Package 'FLicc' is required for fitFLICC().")

  vb <- .vonbFromStock(object)
  mat50 <- mean(.mat50FromStock(object), na.rm=TRUE)
  lor <- .lorFromStock(object)
  wl <- .wlFromStock(object, linf=vb[["linf"]], k=vb[["k"]], t0=vb[["t0"]])

  pars <- FLCore::FLPar(c(linf=vb[["linf"]], k=vb[["k"]], t0=vb[["t0"]], mat50=mat50, lor))
  pars <- rbind(pars, .lorToGis(pars))

  if (is.null(lengthData)) {
    if (!requireNamespace("icesdata", quietly=TRUE))
      stop("Provide 'lengthData' or install package 'icesdata' for lenSamples().")
    lengthData <- icesdata::lenSamples(object, pars, nLengths)
  }

  llb <- as.numeric(dimnames(lengthData)[[1]])
  l50 <- .vbFn(age=mat50, linf=vb[["linf"]], k=vb[["k"]], t0=vb[["t0"]])

  dl <- list(
    model_name=modelName,
    LLB=llb,
    fq=list(swt=lengthData[, 1]),
    Linf=c(as.numeric(pars["linf"]), 0.01),
    sel_fun=selFun,
    Catch=c(catch),
    gear_names="gear1",
    Mk=mk,
    ref_length=refLength,
    a=c(unname(wl["a"])),
    b=c(unname(wl["b"])),
    L50=c(l50),
    GL=c(l50 * 0.9)
  )

  fit <- FLicc::fiticc(dl, compile=compile)

  list(parameters=pars, input=dl, fit=fit)
})

#' @rdname fitFLICC
#' @exportMethod fitFLICC
setMethod("fitFLICC", signature(object="FLStocks"),
function(object, ...) {
  lapply(object, function(stk) fitFLICC(stk, ...))
})

#' FishLife Parameter Interface by Stock
#'
#' Replicates the prototype `ParFishlife <- lapply(split(spp, spp$stock), ...)`
#' workflow in a reusable function.
#'
#' @param spp A data.frame with at least stock, genus, species columns.
#' @param stockCol Name of stock column.
#' @param genusCol Name of genus column.
#' @param speciesCol Name of species column.
#'
#' @return A named list of \code{FLPar} (`linf`, `k`, `t0`, `l50`, `m`) by stock.
#'   If no FishLife match is found for a stock, values are returned as `NA`.
#' @export
fishLifeParsByStock <- function(spp,
                                stockCol="stock",
                                genusCol="Genus",
                                speciesCol="Species") {
  if (!requireNamespace("FishLife", quietly=TRUE))
    stop("Package 'FishLife' is required for fishLifeParsByStock().")
  if (!all(c(stockCol, genusCol, speciesCol) %in% names(spp)))
    stop("spp must contain columns: ", paste(c(stockCol, genusCol, speciesCol), collapse=", "))

  splitSpp <- split(spp, spp[[stockCol]])
  lapply(splitSpp, function(tt) {
    taxa <- suppressMessages(tryIt(
      FishLife::Search_species(
        Genus=as.character(tt[[genusCol]][1]),
        Species=as.character(tt[[speciesCol]][1])
      )$match_taxonomy
    ))
    if (is.null(taxa)) {
      taxa <- suppressMessages(tryIt(
        FishLife::Search_species(
          Genus=as.character(tt[[genusCol]][1])
        )$match_taxonomy
      ))
    }
    if (is.null(taxa))
      return(FLCore::FLPar(c(linf=NA_real_, k=NA_real_, t0=NA_real_, l50=NA_real_, m=NA_real_)))

    pred <- suppressMessages(tryIt(FishLife::Plot_taxa(taxa, Plot=FALSE)))
    if (is.null(pred))
      pred <- suppressMessages(tryIt(FishLife::PlotTaxa(taxa, Plot=FALSE)))
    if (is.null(pred))
      return(FLCore::FLPar(c(linf=NA_real_, k=NA_real_, t0=NA_real_, l50=NA_real_, m=NA_real_)))

    meanPred <- pred[[1]]$Mean_pred
    if (is.null(meanPred))
      return(FLCore::FLPar(c(linf=NA_real_, k=NA_real_, t0=NA_real_, l50=NA_real_, m=NA_real_)))

    # FishLife predictions are on log scale for core life-history metrics.
    mexp <- exp(meanPred)
    FLCore::FLPar(c(
      linf=unname(mexp["Loo"]),
      k=unname(mexp["K"]),
      t0=NA_real_,
      l50=unname(mexp["Lm"]),
      m=unname(mexp["M"])
    ))
  })
}

#' @rdname fishLifeParsByStock
#' @export
fishlifeParsByStock <- fishLifeParsByStock
