# =============================================================================
# Rebuild Methods for biodyn and FLPar Objects
# =============================================================================
# Note: Generic definitions are in generic.R

validPars <- function(object) {
  # Check if required parameters exist
  parNms=dimnames(object)$params
  if (!all(c("r",  "p") %in% parNms) & !all(c("msy","bmsy") %in% parNms)){
    return(NULL)
    warning("'r and p' or 'msy and bmsy', not in object")}
  
  if ("k"%in%parNms & !"virgin"%in%parNms)
    dimnames(object)$params[dimnames(object)$params=="k"]="virgin"
  if (!any(c("virgin") %in% parNms))
    object=rbind(object,FLPar(virgin=1))
  if (any(c(object["virgin"]) <= 0, na.rm = TRUE)) 
    stop("virgin parameter must be positive")
  
  if (!all(c("r", "p") %in% parNms) & all(c("msy","bmsy") %in% parNms)){
    if (any(c(object["msy"]) <= 0, na.rm = TRUE)) 
      stop("msy parameter must be positive")
    if (any(c(object["bmsy"]) <= 0, na.rm = TRUE)) 
      stop("bmsy parameter must be positive")
    
    if (!("p"%in%dimnames(object)$params))
       object=rbind(object,p(object))  
    if (!("r"%in%dimnames(object)$params))
      object=rbind(object,r(object))}
  
  if (any(c(object["r"]) <= 0, na.rm = TRUE)) 
    stop("r parameter must be positive")
  
  if (!("bmsy"%in%dimnames(object)$params)) 
    object=rbind(object,bmsy(object))
  if (!("msy"%in%dimnames(object)$params)) 
    object=rbind(object,msy(object))
  
  return(object)}

# rebuildTime generic is defined in generic.R

rebuildTimeFn <- function(object, nYrs = 50, minVal = 0.05) {
  
  y = biodyn(object, nyrs = nYrs)
  
  y@stock[]=refpts(y)["bmsy"]*minVal
  
  z=fwd(y,catch=catch(y)[,-1]%=%0)
  
  df=as.data.frame(stock(z)%/%refpts(y)["bmsy"],drop=T)[-1,]
  df=subset(transmute(df,initial=data,T_recover=subset(df,data>1)[1,"year"]-year),initial<=1)
  
  rbind(df[order(df$initial),],data.frame(initial=1,T_recover=0))}


#' Calculate Rebuild Time using FLPar parameters
#'
#' @description A simplified and faster version of rebuildTime that works directly with FLPar objects.
#' This function calculates rebuilding time using Pella-Tomlinson dynamics parameters.
#'
#' @param object An FLPar object containing parameters (r, k, p, B0)
#' @param nyrs Number of years for projection (default = 50)
#' @param minVal Minimum depletion value as proportion of BMSY (default = 0.05)
#' @param ... Additional arguments
#' @return A data frame with columns:
#'   \item{initial}{Initial depletion level relative to BMSY}
#'   \item{T_recover}{Time to recover to BMSY}
#' @export
#' @examples
#' # Create FLPar with Pella-Tomlinson parameters
#' params=FLPar(r = 0.5, k = 1000, p = 1, B0 = 1000)
#' rebuild_time=rebuildTime2(params)
#' 
#' # With custom parameters
#' rebuild_time=rebuildTime2(params, nyrs = 100, minVal = 0.1)
# rebuildTime2 and rebuildTime3 generics are defined in generic.R

#' @rdname rebuildTime2
#' @export
setMethod("rebuildTime2", signature(object = "FLPar"),
  function(object, nyrs = 50, ftar = 0.0, timing = 0,  minVal = 0.05, ...) {
    
    # Validate other arguments
    if (!is.numeric(nyrs) || nyrs <= 0)
      stop("nyrs must be a positive integ/er")
    if (!is.numeric(ftar) || ftar < 0) 
      stop("ftar must be non-negative")
    if (!is.numeric(timing) || timing < 0 || timing>1) 
      stop("timimg must be [0,1]")

    object=validPars(object)
    
    dimnames(object)$params[dimnames(object)$params=="virgin"]="k"
    
    # Call the existing rebuildTimeFn function (avoiding recursion)
    y=biodyn(object, nyrs = nyrs)
    
    y@stock[]=refpts(y)["bmsy"] * minVal
    
    z=mpb::fwd(y, catch = catch(y)[,-1] %=% 0)
    
    df=as.data.frame(stock(z) %/% refpts(y)["bmsy"], drop = TRUE)[-1,]
    df=subset(transmute(df, initial = data, T_recover = subset(df, data > 1)[1,"year"] - year), initial <= 1)
    df=rbind(df[order(df$initial),],data.frame(initial=1,T_recover=0))
  
    return(df)})

#' Calculate Rebuild Time using FLPar parameters (Alternative Implementation)
#'
#' @description An alternative implementation of rebuildTime that works directly with FLPar objects.
#' This version uses a different forward projection approach with explicit surplus production calculations.
#'
#' @param object An FLPar object containing parameters (r, p, bmsy, fmsy, virgin)
#' @param nyrs Number of years for projection (default = 50)
#' @param ftar Target fishing mortality multiplier (default = 0.0 for no fishing)
#' @param timing Timing parameter (default = 0)
#' @param ... Additional arguments
#' @return A data frame with columns:
#'   \item{initial}{Initial depletion level relative to BMSY}
#'   \item{T_recover}{Time to recover to BMSY}
#' @export
#' @examples
#' # Create FLPar with required parameters
#' params=FLPar(r = 0.5, p = 1, bmsy = 500, fmsy = 0.2, virgin = 1000)
#' rebuild_time=rebuildTime3(params)
#' 
#' # With fishing mortality
#' rebuild_time=rebuildTime3(params, ftar = 0.5)
setMethod("rebuildTime3", signature(object = "FLPar"),
  function(object, nyrs = 50, ftar = 0.0, timing = 0, ...) {
    
    # Validate other arguments
    if (!is.numeric(nyrs) || nyrs <= 0)
      stop("nyrs must be a positive integ/er")
    if (!is.numeric(ftar) || ftar < 0) 
      stop("ftar must be non-negative")
    if (!is.numeric(timing) || timing < 0 || timing>1) 
      stop("timimg must be [0,1]")
    
    object=validPars(object)
    
    # Call the forward projection function
    return(fwdPT(object[c("r","p","virgin","bmsy","msy")], nyrs=nyrs, ftar=ftar, timing=timing, ...))})

#' Forward Projection Function for Rebuild Time Calculation
#'
#' @description Internal function that performs forward projection using Pella-Tomlinson dynamics.
#' This is the core calculation engine for rebuildTime3.
#'
#' @param object FLPar object with parameters
#' @param nyrs Number of years for projection
#' @param ftar Target fishing mortality multiplier
#' @param timing Timing parameter
#' @param ... Additional arguments
#' @return Data frame with rebuild time results
#' @keywords internal
fwdPT=function(object, nyrs = 50, ftar = 0.0, timing = 0) {
  
  nits=dim(object)[2]
  stock  =FLQuant(0.05,dimnames=list(year=seq(nyrs),iter=seq(nits)))%*%object["bmsy"]
  harvest=FLQuant(1,   dimnames=list(year=seq(nyrs),iter=seq(nits)))%*%object["msy"]%/%object["bmsy"]*ftar
  catch  =FLQuant(0,   dimnames=list(year=seq(nyrs),iter=seq(nits)))
  sp     =FLQuant(0,   dimnames=list(year=seq(nyrs),iter=seq(nits))) 
  
  for (y in dimnames(stock)$year[-dim(stock)[2]]) {
    
    sp[, y]=(object["r"] %/% object["p"]) %*% stock[, y] %*% (1 - exp(log(stock[, y] / object["virgin"]) %*% object["p"]))
    catch[, y]=stock[, y] * harvest[, y]
    stock[, ac(an(y) + 1)]=stock[, y] - catch[, y] + sp[, y]}
  
  stock=stock %/% object["bmsy"]

  df=as.data.frame(stock, drop = TRUE)
  df=subset(transmute(df, initial = data, T_recover = subset(df, data > 1)[1, "year"] - year), initial <= 1)
  
  rbind(df[order(df$initial), ], data.frame(initial = 1, T_recover = 0))
}

#' Calculate Rebuild Trajectories
#'
#' @description Projects stock rebuilding trajectories from different initial depletion levels
#'
#' @param object A biodyn object
#' @param target Target biomass (default = BMSY)
#' @param nInitial Number of initial depletion levels (default = 100)
#' @param growthRate Growth rate for depletion sequence (default = 0.3)
#' @param minVal Minimum depletion value (default = 1e-6)
#' @param maxVal Maximum depletion value (default = 1)
#' @param nx Number of interpolation points (default = 101)
#' @return A data frame with columns:
#'   \item{year}{Projection year}
#'   \item{initial}{Initial depletion level relative to BMSY}
#' @export
#' @examples
#' bd=biodyn(FLPar(r=0.5, k=1000, p=1))
#' rebuild_data=rebuild(bd)
# rebuildTime generic with extended signature is defined in generic.R
# This signature variation is handled via method dispatch

setMethod("rebuildTime", signature(object="biodyn"), 
          function(object, target=refpts(object)["bmsy"], nInitial=100, 
                   growthRate=0.3, minVal=1e-6, maxVal=1, nx=101) {
            if (!is(object, "biodyn"))
              stop("object must be a biodyn object")
            
            if (!is.numeric(target) || length(target) != 1)
              stop("target must be a single numeric value")
            
            if (!all(sapply(list(nInitial, nx), function(x) is.numeric(x) && x > 0)))
              stop("nInitial, and nx must be positive integers")
            
            if (!all(sapply(list(growthRate, minVal, maxVal), is.numeric)))
              stop("growthRate, minVal, and maxVal must be numeric")
            
            minVal=1e-6
            if (minVal >= maxVal)
              stop("minVal must be less than maxVal")
            
            bmsy=c(refpts(object)["bmsy"])
            
            rtn=propagate(object, nInitial)
            
            # Create stock projection
            target_seq=c(target)*seq(minVal^growthRate, maxVal^growthRate, length.out=nInitial)^(1/growthRate)
            rtn@stock=FLQuant(rep(target_seq, each=dim(rtn)[2]), dimnames=dimnames(stock(rtn)))
            rtn=fwd(rtn, catch=catch(rtn)[,-1]%=%0.0)
            
            # Transform data
            dat=as.data.frame(stock(rtn), drop=TRUE)
            dat$initial=c(stock(rtn)[,1])[an(dat$iter)]
            dat=dat[,-2]
            
            # Interpolate results
            dat=as.data.frame(with(dat, akima::interp(initial, 
                                                      data, year, yo=bmsy, duplicate="mean", nx=nx, jitter=1e-6)))[,c(3,1)]
            names(dat)=c("year", "initial")
            
            transform(dat,initial=initial/bmsy)})

#' @rdname rebuild
#' @export
setMethod("rebuild", signature(object="numeric"),
  function(object, ...) {
    args = list(...)
    if (is.null(args$p)) stop("Argument 'p' must be supplied in ...")
    p = args$p
    k = if (!is.null(args$k)) args$k else 1e3
    b0 = if (!is.null(args$b0)) args$b0 else 1
    nyrs = if (!is.null(args$nyrs)) args$nyrs else 50
    niters = if (!is.null(args$niters)) args$niters else 101
    r = object
    object = biodyn(FLPar(r = r, p = p, k = k, B0 = b0))
    shape = c(refpts(object)["bmsy"] %/% params(object)["k"])
    object = window(object, end = nyrs)
    object@stock[] = refpts(object)["bmsy"]
    object@catch[] = 0
    target = c(refpts(object)["bmsy", 1])
    object = propagate(object, niters)
    object@stock = object@stock %*% FLQuant(rep(seq(0, 1, length.out = niters), each = nyrs), dimnames = dimnames(stock(object)))
    object = fwd(object, harvest = stock(object)[, -1] %=% 0)
    dat = cbind(target = target, as.data.frame(stock(object), drop = TRUE))
    
    # Use data.table for fast grouping and min search
    dt = data.table::as.data.table(dat)
    dt[, diff2 := (data - target)^2]
    dt_min = dt[, .SD[which.min(diff2)], by = iter]
    
    # Add initial
    dt_min[, initial := c(stock(object)[,1,,,,iter]) / target]
    dt_min = dt_min[order(initial), .(year, initial)]
    dt_min = dt_min[-1, ]
    out = cbind(shape = shape, dt_min[year < nyrs, ])
    
    as.data.frame(out)})

