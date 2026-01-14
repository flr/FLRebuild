# =============================================================================
# Stock-Recruitment Relationship Functions
# =============================================================================
# Note: gsl package should be in package dependencies

# =============================================================================
# S4 Generics
# =============================================================================

#' Generic for inverse stock-recruitment relationships
#' 
#' @param params FLPar object with SRR parameters or FLBRP/FLSR object
#' @param rec Recruitment values (numeric, FLPar, or FLQuant)
#' @param ... Additional arguments
#' @return SSB values corresponding to the given recruitment
#' @export
setGeneric("invSRR", function(params, rec, ...) standardGeneric("invSRR"))

#' Create reference point FLPar object
#' 
#' @param object An object (typically FLBRP)
#' @param ref Reference point name (character)
#' @param value Value for the reference point (numeric, default NA)
#' @param quant Quant to set (default is ref)
#' @param ... Additional arguments
#' @return FLPar object with reference point structure
#' @export
setGeneric("refCreate", function(object, ref, value=NA, quant=ref, ...) standardGeneric("refCreate"))

#' Maximum recruitment for a stock-recruitment model
#' 
#' @param object An object (typically FLBRP)
#' @param ratio Optional numeric ratio (default missing)
#' @param ... Additional arguments
#' @return Maximum recruitment (numeric or FLPar)
#' @export
setGeneric("rmax", function(object, ratio, ...) standardGeneric("rmax"))

#' Recruitment at MSY for a stock-recruitment model
#' 
#' @param object An object (typically FLBRP)
#' @param ratio Optional numeric ratio (default 1.0)
#' @param ... Additional arguments
#' @return Recruitment at MSY (numeric or FLPar)
#' @export
setGeneric("rmsy", function(object, ratio=1.0, ...) standardGeneric("rmsy"))


#' Recruitment at Virgin for a stock-recruitment model
#' 
#' @param object An object (typically FLBRP)
#' @param ratio Optional numeric ratio (default 1.0)
#' @param ... Additional arguments
#' @return Recruitment at MSY (numeric or FLPar)
#' @export
setGeneric("rvirgin", function(object, ratio=1.0, ...) standardGeneric("rvirgin"))

# =============================================================================
# S4 Methods
# =============================================================================

#' @rdname invSRR
#' @export
setMethod("invSRR", signature(params="FLBRP", rec="FLQuant"),
  function(params, rec, ...) {
    switch(SRModelName(model(params)),
           bevholt=bevholtInv(params(params), rec),
           ricker =rickerInv( params(params), rec))
  })

#' @rdname invSRR
#' @export
setMethod("invSRR", signature(params="FLBRP", rec="FLPar"),
  function(params, rec, ...) {
    switch(SRModelName(model(params)),
           bevholt=bevholtInv(params(params), rec),
           ricker=rickerInv(params(params), rec))
  })

#' @rdname invSRR
#' @export
setMethod("invSRR", signature(params="FLSR", rec="FLQuant"),
  function(params, rec, ...) {
    switch(SRModelName(model(params)),
           bevholt=bevholtInv(params(params), rec),
           ricker=rickerInv(params(params), rec))
  })

#' @rdname refCreate
#' @details
#' For an FLBRP object, returns a new FLPar reference point object with the specified values. The input FLBRP object is not modified.
#' For an FLPar object, returns a new FLPar reference point object mapping parameter names to quant columns. The input FLPar object is not modified.
#' @return An FLPar object with reference point structure.
#' @export
setMethod("refCreate", signature(object="FLBRP"), 
  function(object, ref, value=NA, quant=ref, ...) {
    # Validate quant
    valid_quants=c("harvest", "yield", "rec", "ssb", "biomass", "revenue", "cost", "profit")
    if (any(!quant %in% valid_quants)) {
      stop("Invalid quant(s) provided. Must be one of: ", paste(valid_quants, collapse=", "))
    }
    # Determine iter dimension
    niter=if (inherits(value, "FLPar") || inherits(value, "FLQuant")) {
      if (!is.null(dim(value)["iter"])) dim(value)["iter"] else 1
    } else 1
    # Build dimnames
    dn=list(refpt=ref, quant=valid_quants, iter=seq_len(niter))
    refpar=FLPar(NA, dimnames=dn)
    # Assign value(s) to specified quant(s)
    if (length(quant) == 1 && length(value) > 1) {
      # If a single quant but multiple values, assign across iter
      refpar[ref, quant, ]=value
    } else {
      # Otherwise assign value(s) to quant(s)
      refpar[ref, quant]=value
    }
    return(refpar)
  })

setMethod("refCreate", signature(object="FLPar"), 
  function(object, ...) {
    dmns=dimnames(object)
    valid_quants=c("harvest", "yield", "rec", "ssb", "biomass", "revenue", "cost", "profit")
    niter=if (!is.null(dmns$iter)) length(dmns$iter) else 1
    dn=list(refpt=dmns$params, quant=valid_quants, iter=dmns$iter)
    refpar=FLPar(NA, dimnames=dn)
    # Map parameter names to quant columns
    if (any(startsWith(dmns$params, "b")))
      refpar[startsWith(dmns$params, "b"), "ssb", ]=object[startsWith(dmns$params, "b"), ]
    if (any(startsWith(dmns$params, "f")))
      refpar[startsWith(dmns$params, "f"), "harvest", ]=object[startsWith(dmns$params, "f"), ]
    return(refpar)})

setMethod("refCreate", signature(object="character"), 
          function(object, ref="", value=NA, quant="ssb", ...) {
          
            dmns=dimnames(refpts(FLBRP()))
            dmns$refpt=object
            
            rtn=FLPar(NA,dimnames=dmns)
            rtn[,quant]=value
            
            rtn})
            
#' @rdname rmax
#' @export
setMethod("rmax", signature(object="FLBRP", ratio="missing"), 
  function(object, ...) {
    rtn=switch(SRModelName(model(object)),
           bevholt = params(object)["a"],
           ricker = 1 / params(object)["b"],
           segreg = params(object)["a"] %/% params(object)["a"])
    rtn})

#' @rdname rmax
#' @export
setMethod("rmax", signature(object="FLBRP", ratio="numeric"), 
  function(object, ratio, ...) {
    if (ratio <= 0 || ratio > 1) return(refCreate(object, "rmax", NA, "ssb"))
    
    # Create new reference point FLPar
    rmaxPar=refCreate(object, paste("rmax", ratio*100, sep="."), invSRR(object, rmax(object) * ratio), "ssb")
    # Add to refpts slot
    refpts(object)=rmaxPar
    rtn=computeRefpts(object)
    
    rtn})

#' @rdname rmsy
#' @export
setMethod("rmsy", signature(object="FLBRP", ratio="missing"), 
  function(object, ratio=1.0, ...) {
    rec_val=FLPar(refpts(object)["msy", "rec", drop=TRUE]) * ratio
    rmsyPar=refCreate(object, "rmsy", invSRR(object, rec_val), "rec")
    refpts(object)=rbind(refpts(object), rmsyPar)
    computeRefpts(object)
  })

#' @rdname rmsy
#' @export
setMethod("rmsy", signature(object="FLBRP", ratio="numeric"), 
  function(object, ratio, ...) {
    if (ratio <= 0 || ratio > 1) return(refCreate(object, "rmsy", NA, "ssb"))
    
    refpts(object)=refCreate(object, "msy", NA, "ssb")
    refpts(object)=computeRefpts(object)
    rec_val=FLPar(refpts(object)["msy", "rec", drop=TRUE]) * ratio
    
    refpts(object)=refCreate(object, paste("rmsy", ratio*100, sep="."), invSRR(object, rec_val), "rec")
    rtn=computeRefpts(object)
    
    rtn})

#' @rdname rvirgin
#' @export
setMethod("rvirgin", signature(object="FLBRP", ratio="missing"), 
  function(object, ratio=1.0, ...) {
    rec_val=FLPar(refpts(object)["virgin", "rec", drop=TRUE]) * ratio
    rvirginPar=refCreate(object, "rvirgin", invSRR(object, rec_val), "rec")
    refpts(object)=rbind(refpts(object), rvirginPar)
    computeRefpts(object)
  })

#' @rdname rvirgin
#' @export
setMethod("rvirgin", signature(object="FLBRP", ratio="numeric"), 
  function(object, ratio, ...) {
    if (ratio <= 0 || ratio > 1) return(refCreate(object, "rvirgin", NA, "ssb"))
    
    rec_val=FLPar(refpts(object)["virgin", "rec", drop=TRUE]) * ratio
    refpts(object)=refCreate(object, paste("rvirgin", ratio*100, sep="."), invSRR(object, rec_val), "rec")
    rtn=computeRefpts(object)
    
    rtn})

# =============================================================================
# Helper Functions (not S4)
# =============================================================================

#' Inverse Beverton-Holt stock-recruitment relationship
#' 
#' @param params FLPar with parameters 'a' and 'b'
#' @param rec Recruitment values
#' @return SSB values
#' @keywords internal
bevholtInv=function(params, rec) {
  return((rec %*% params["b"]) %/% (params["a"] %-% rec))
}

#' Inverse Ricker stock-recruitment relationship
#' 
#' @param params FLPar with parameters 'a' and 'b'
#' @param rec Recruitment values
#' @return SSB values
#' @keywords internal
rickerInv=function(params, rec) {
  arg=-(rec %*% params["b"]) %/% params["a"]
  arg[arg > 0]=NA
  arg[arg <= 0]=sapply(arg[arg <= 0], gsl::lambert_Wm1)
  SSB=-1 %/% (params["b"] %*% arg)
  return(SSB)
}

#' Maximum recruitment SSB for Ricker model
#' 
#' @param b Ricker parameter b
#' @return SSB at maximum recruitment
#' @keywords internal
rickerMaxRec=function(b) return(1/b)


