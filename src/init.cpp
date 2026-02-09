// Combined initialization for FLRebuild
// This file combines TMB and Rcpp initializations to avoid duplicate R_init_FLRebuild definitions

#include <R.h>
#include <R_ext/Rdynload.h>

// Forward declaration of TMB init function (defined by TMB via TMB_LIB_INIT macro)
extern "C" void R_init_FLRebuild_TMB(DllInfo *dll);

// Forward declaration of Rcpp function (from RcppExports.cpp)
extern "C" SEXP _FLRebuild_loglAR1(SEXP obsSEXP, SEXP hatSEXP, SEXP rhoSEXP);

// Rcpp function registration table
static const R_CallMethodDef CallEntries[] = {
    {"_FLRebuild_loglAR1", (DL_FUNC) &_FLRebuild_loglAR1, 3},
    {NULL, NULL, 0}
};

// Main init function that combines both TMB and Rcpp registrations
extern "C" void R_init_FLRebuild(DllInfo *dll) {
    // Register TMB functions (this will register all TMB-generated functions)
    R_init_FLRebuild_TMB(dll);
    
    // Register Rcpp functions
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
