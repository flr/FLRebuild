# Package initialization
# This file is loaded last (due to alphabetical ordering) and handles DLL loading
# Following the standard TMB pattern used in FLCandy

.onLoad <- function(libname, pkgname) {
  # Load the TMB DLL when the package is loaded
  # Standard TMB approach: use dyn.load(TMB::dynlib("DLLname"))
  # Try FLCandy first (legacy DLL that previously worked), then FLRebuild
  
  # Try FLCandy first (legacy)
  tryCatch({
    dyn.load(TMB::dynlib("FLCandy"))
  }, error = function(e) {
    # If FLCandy not found, try FLRebuild
    tryCatch({
      dyn.load(TMB::dynlib("FLRebuild"))
    }, error = function(e2) {
      # Silent failure - will be caught by .ensureTMBdll() when MakeADFun is called
    })
  })
}

.onUnload <- function(libpath) {
  # Unload the DLLs when package is unloaded
  tryCatch({
    dll <- getLoadedDLLs()[["FLCandy"]]
    if (!is.null(dll)) library.dynam.unload("FLCandy", libpath)
  }, error = function(e) NULL)
  tryCatch({
    dll <- getLoadedDLLs()[["FLRebuild"]]
    if (!is.null(dll)) library.dynam.unload("FLRebuild", libpath)
  }, error = function(e) NULL)
}
