# Package initialization
# This file is loaded last (due to alphabetical ordering) and handles DLL loading

.onLoad <- function(libname, pkgname) {
  # Load the TMB DLL when the package is loaded
  # This ensures the DLL is available when MakeADFun is called
  
  # Try to get DLL path using TMB::dynlib (works for installed packages)
  dll_file <- tryCatch({
    TMB::dynlib("FLRebuild")
  }, error = function(e) {
    # Fallback: construct path manually
    system.file("libs", paste0(pkgname, .Platform$dynlib.ext), package = pkgname)
  })
  
  # Also try src directory for development installations
  if (!file.exists(dll_file)) {
    dll_file <- system.file("src", paste0(pkgname, .Platform$dynlib.ext), package = pkgname)
  }
  
  if (file.exists(dll_file)) {
    dyn.load(dll_file, local = FALSE, now = TRUE)
  } else {
    warning("FLRebuild DLL not found at expected locations. ",
            "The package may need to be recompiled and reinstalled.")
  }
}

.onUnload <- function(libpath) {
  # Unload the DLL when package is unloaded
  library.dynam.unload("FLRebuild", libpath)
}
