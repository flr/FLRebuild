// Combined initialization for FLRebuild package
// Handles both TMB and Rcpp initialization

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Forward declarations
extern void R_init_FLRebuild_TMB(DllInfo *dll);
extern void R_init_FLRebuild_Rcpp(DllInfo *dll);

// Combined initialization function
extern "C" void R_init_FLRebuild(DllInfo *dll) {
  // Initialize TMB components
  R_init_FLRebuild_TMB(dll);
  
  // Initialize Rcpp components  
  R_init_FLRebuild_Rcpp(dll);
}
