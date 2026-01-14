# =============================================================================
# Analysis and Results Functions
# =============================================================================
# Note: Generic definition is in 00-generic.R

#' @rdname results
#' @export
setMethod("results", signature(object = "FLStock", eq = "FLBRP"),
  function(object, eq, ftar = 0.5) {
    
    # Calculate relative SSB
    ssb_rel <- ssb(object) %/% refpts(eq)["msy", "ssb"]
    
    # Calculate relative recruits (age 0 stock numbers / virgin recruits)
    rec_rel <- stock.n(object)[1, ]
    dimnames(rec_rel)$age <- "all"
    rec_rel <- rec_rel %/% refpts(eq)["virgin", "rec"]
    
    # Calculate relative yield
    yield_rel <- catch(object) %/% refpts(eq)["msy", "yield"]
    
    # Calculate relative F
    f_rel <- fbar(object) %/% refpts(eq)["msy", "harvest"]
    
    # Calculate relative forage (need to set fbar on eq for MSY forage)
    eq_forage <- eq
    fbar(eq_forage) <- fbar(eq_forage)[, 1]
    fbar(eq_forage)[] <- refpts(eq_forage)["msy", "harvest"]
    forage_rel <- forage(object) %/% forage(eq_forage)
    
    # Calculate ABI
    abi_vals <- abi(object, eq)
    
    # Calculate rebuild metric
    rebuild_vals <- pBuild(ssb(object), eq, 
                           ftar = ftar * refpts(eq)["msy", "harvest"])
    names(dimnames(rebuild_vals))[1] <- "age"
    
    # Combine into FLQuants
    rtn <- FLQuants(
      SSB      = ssb_rel,
      Recruits = rec_rel,
      Yield    = yield_rel,
      F        = f_rel,
      Forage   = forage_rel,
      ABI      = abi_vals,
      Rebuild  = rebuild_vals
    )
    
    return(rtn)
  })

#' @rdname results
#' @export
setMethod("results", signature(object = "list", eq = "FLBRP"),
  function(object, eq, ftar = 0.5) {
    # Handle list of FLStock objects (e.g., from MSE results)
    # Apply results method to each element
    lapply(object, function(x) {
      if (inherits(x, "FLStock")) {
        results(x, eq, ftar = ftar)
      } else {
        x
      }
    })
  })
