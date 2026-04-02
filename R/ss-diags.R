#' Run Stock Synthesis diagnostics
#'
#' Run a standard set of diagnostics from Stock Synthesis model output,
#' optionally generating figures and saving the loaded report.
#'
#' @param object Directory containing Stock Synthesis output files.
#' @param name Character label used for output filenames.
#' @param plot Logical; if `TRUE`, create diagnostic plots.
#'
#' @return A list with components `model` (the `SS_output` object) and `mvln`
#'   (the result from `SSdeltaMVLN`).
#' @export
setMethod("ssDiagnostics", signature(object = "character"), function(object, name, plot = TRUE) {
    # Create diagnostics plot directory if plotting enabled
    if(plot) dir.create("Plotdiags", showWarnings = FALSE)
    
    # Load model
    ss3rep <- r4ss::SS_output(dir = object, covar = TRUE)
    save(ss3rep, file = paste0(name, ".Rdata"))
    
    if(plot) {
      # Basic R4SS plots
      r4ss::SS_plots(ss3rep, dir = getwd())
      
      # Data setup plots
      r4ss::sspar()
      r4ss::SSplotData(ss3rep, subplots = 2)
      dev.print(jpeg, paste0("Plotdiags/DataSetup_", name, ".jpg"), 
                width = 8, height = 6, res = 300, units = "in")
      
      # Residual diagnostics
      r4ss::sspar(mfrow = c(2,4), plot.cex = 0.8)
      r4ss::SSplotRunstest(ss3rep, subplots = "cpue", add = TRUE, legendcex = 0.6, mixing = "two.sided")
      r4ss::SSplotRunstest(ss3rep, subplots = "age",  add = TRUE, legendcex = 0.6, mixing = "less")
      dev.print(jpeg, paste0("Plotdiags/RunsTestResiduals_", name, ".jpg"), 
                width = 8, height = 7, res = 300, units = "in")
      
      # Joint residuals
      r4ss::sspar(mfrow = c(1,2), plot.cex = 0.8)
      ss3diags::SSplotJABBAres(ss3rep, subplots = "cpue", add = TRUE, col = r4ss::sscol(3)[c(1,3,2)])
      ss3diags::SSplotJABBAres(ss3rep, subplots = "age", add = TRUE, col = r4ss::sscol(3)[c(1,3,2)])
      dev.print(jpeg, paste0("Plotdiags/JointResiduals_", name, ".jpg"), 
                width = 8, height = 3.5, res = 300, units = "in")
      
      # MVLN uncertainty
      r4ss::sspar(mfrow = c(1,1), plot.cex = 0.9)
      mvn <- r4ss::SSdeltaMVLN(ss3rep, plot = TRUE, Fref = c("Btgt"), catch.type = c("Exp"))
      kbproj <- data.frame(mvn$kb)
      r4ss::SSplotKobe(kbproj, fill = TRUE, joint = FALSE, posterior = "kernel",
                 ylab = expression(F/F[trg]), xlab = expression(SSB/SSB[trg]))
      dev.print(jpeg, paste0("Plotdiags/Kobe_", name, ".jpg"), 
                width = 6.5, height = 6.5, res = 300, units = "in")
    } else {
      mvn <- r4ss::SSdeltaMVLN(ss3rep, plot = FALSE, Fref = c("Btgt"), catch.type = c("Exp"))
    }
    
    return(list(model = ss3rep, mvln = mvn))
})

#' Run Stock Synthesis retrospective diagnostics
#'
#' Summarize retrospective model runs and optionally create standard
#' retrospective and hindcast cross-validation plots.
#'
#' @param object A list of retrospective model outputs.
#' @param name Character label used for output filenames.
#' @param plot Logical; if `TRUE`, create diagnostic plots.
#'
#' @return A list with components `retroSummary` and `hccomps`.
#' @export
setMethod("ssRetrospective", signature(object = "list"), function(object, name, plot = TRUE) {
  retroSummary <- r4ss::SSsummarize(object)
  hccomps <- r4ss::SSretroComps(object)
  
  if(plot) {
    # Retrospective plots
    r4ss::sspar(mfrow = c(2,2), plot.cex = 0.9)
    r4ss::SSplotRetro(retroSummary, forecastrho = TRUE, add = TRUE, subplots = "SSB")
    r4ss::SSplotRetro(retroSummary, forecastrho = TRUE, add = TRUE, legend = FALSE)
    r4ss::SSplotRetro(retroSummary, subplots = "F", add = TRUE, legendloc = "left", legendcex = 0.8)
    r4ss::SSplotRetro(retroSummary, subplots = "F", forecastrho = TRUE, add = TRUE, legend = FALSE)
    dev.print(jpeg, paste0("Plotdiags/RetroForecast_", name, ".jpg"), 
              width = 8, height = 9, res = 300, units = "in")
    
    # Hindcast cross-validation
    r4ss::sspar(mfrow = c(1,2), plot.cex = 0.9)
    r4ss::SSplotHCxval(retroSummary, add = TRUE, legendcex = 0.6, Season = 1)
    dev.print(jpeg, paste0("Plotdiags/HCxvalIndex_", name, ".jpg"), 
              width = 8, height = 5, res = 300, units = "in")
  }
  
  return(list(retroSummary = retroSummary, hccomps = hccomps))
})

# For basic diagnostics
#model=ssDiagnostics(path="P:/rfmo/ices/wkbseabass/ss3/north/base",name="test")

# For retrospective analysis (if you have retro models)
#retro_results=ssRetrospective(retros=, name="test")