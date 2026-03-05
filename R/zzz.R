# Package initialization
# This file is loaded last (due to alphabetical ordering)
# Note: DLL loading is handled by useDynLib() in NAMESPACE
# TMB_LIB_INIT in src/FLSRTMB.cpp causes TMB to automatically create R_init_FLRebuild
# which registers TMB symbols like getParameterOrder
# Rcpp routines should be registered via attributes or will be handled separately
