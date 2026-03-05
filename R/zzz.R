# Package initialization
# This file is loaded last (due to alphabetical ordering) and handles DLL loading

.onLoad <- function(libname, pkgname) {
  # Load the TMB DLL when the package is loaded
  # This ensures the DLL is available when MakeADFun is called
  
  # Use library.dynam which is the standard R way to load package DLLs
  # This works for both installed packages and development installations
  dll_file <- system.file("libs", paste0(pkgname, .Platform$dynlib.ext), 
                          package = pkgname, lib.loc = libname)
  
  # For development installations, also check src directory
  if (!file.exists(dll_file)) {
    dll_file <- system.file("src", paste0(pkgname, .Platform$dynlib.ext), 
                            package = pkgname, lib.loc = libname)
  }
  
  if (file.exists(dll_file)) {
    # Use library.dynam which properly registers the DLL
    library.dynam("FLRebuild", pkgname, libname, now = TRUE)
  } else {
    # Try using TMB::dynlib as a fallback (may work if package is installed)
    tryCatch({
      dll_file <- TMB::dynlib("FLRebuild")
      if (file.exists(dll_file)) {
        dyn.load(dll_file, local = FALSE, now = TRUE)
      }
    }, error = function(e) {
      warning("FLRebuild DLL not found. The package may need to be recompiled and reinstalled.")
    })
  }
}

.onUnload <- function(libpath) {
  # Unload the DLL when package is unloaded
  library.dynam.unload("FLRebuild", libpath)
}
