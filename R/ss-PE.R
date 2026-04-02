# The figure compares how well different production‐function formulations reproduce the SS3 dynamics in biomass and SSB space, and how much apparent process error remains.
# 
# ## Surplus production panels
# 
# - The top‐left panel shows empirical surplus production against total biomass, overlaid with the fitted Pella–Tomlinson curve in biomass space using a fixed \(B_{\text{MSY}}/K = 0.2\).  
# - The top‐right panel shows empirical surplus production against SSB, overlaid with:
#   - the equilibrium production curve from SS3 (orange), and  
# - the Pella–Tomlinson fit in SSB space with the same \(B_{\text{MSY}}/K\) (red).  
# - Visually, the SSB PT curve tracks the SS3 equilibrium relationship much more closely than the biomass PT curve tracks the biomass SP cloud, indicating that SSB space provides a more consistent production signal for this stock. 
# 
# 
# ## Process‑error panel
# 
# - The bottom panel shows log one‐step‐ahead process error through time for three models:
#   - blue: PT fitted in biomass space (Biomass PT);  
# - red: PT fitted in SSB space (SSB PT);  
# - black dashed: SSB spline PF, i.e. the SS3‐derived production function used in `curveSS` (reference).  
# - The orange series represents the **reference process error** implied by the SS3 equilibrium spline; it retains low‑frequency structure and large swings, highlighting periods where the assessment dynamics depart from a simple production function.  
# - The red series (SSB PT) is close to the orange series in both level and pattern, showing that the PT model in SSB space captures most of the SS3 production dynamics without absorbing all of the apparent process error.  
# - The blue series (Biomass PT) is systematically closer to zero and smoother, demonstrating that fitting PT in biomass space absorbs more of the structural deviation into the production curve itself, leaving less residual process error and therefore providing weaker diagnostics of mis‑specification and prediction skill. 
# 
# 
# Persistent low‑frequency structure and large swings in the process‑error series usually mean the model is **systematically wrong**, not just noisy, and that prediction skill (especially for medium‑term forecasts) is likely to be poor. [edepot.wur](https://edepot.wur.nl/556274)
# 
# ## 1. What low‑frequency structure implies
# 
# When PE is defined as one‑step‑ahead log residuals:
#   
#   - **Trends or regime‑like shifts** in PE indicate that some biological or fishery process (productivity, selectivity, catch reporting, etc.) is changing in a way the model does not represent.  
# - **Strong autocorrelation** (runs of positive or negative PE) shows the model repeatedly under‑ or over‑predicts, again pointing to structural mis‑specification rather than random year‑to‑year variability. [repository.library.noaa](https://repository.library.noaa.gov/view/noaa/52235/noaa_52235_DS1.pdf)
# 
# For your PF‑based diagnostics, that means:
#   
#   - The chosen production function (SS spline or PT) cannot explain the observed biomass/SSB trajectory with stationary noise; unmodelled dynamics are being pushed into PE.
# 
# ## 2. Consequences for prediction skill and advice
# 
# - **Reduced prediction skill**: Hindcasts that rely on a mis‑specified model tend to have biased medium‑term projections, often underestimating risk of decline or overestimating rebuilding speed. [bibbase](https://bibbase.org/network/publication/kell-kirnoto-kitakado-evaluationofthepredictionskillofstockassessmentusinghindcasting-2016)
# - **Biased reference points and status**: If process error compensates for wrong productivity or selectivity assumptions, MSY- and \(B_{\text{MSY}}\)-related quantities can be biased while still giving seemingly reasonable fits. [edepot.wur](https://edepot.wur.nl/556274)
# - **Over‑confident uncertainty**: Standard errors and confidence intervals assume independent, zero‑mean residuals; low‑frequency structure violates this and makes uncertainty bands too narrow, overstating confidence in advice. [spo.nmfs.noaa](https://spo.nmfs.noaa.gov/sites/default/files/TMSPO240.pdf)
# 
# ## 3. How to use this diagnostically
# 
# With your setup:
#   
#   - Treat the SS‑based PF PE (`tseries$pe` / `pe_from_pf`) as the **primary mis‑specification signal**.  
# - Use the similarity or difference between that PE and PT‑based PE (biomass vs SSB) to:
#   - identify which state variable / PF form better preserves the structural signal, and  
# - guide alternative model hypotheses (time‑varying productivity, non‑PT shape, changes in selectivity, etc.). [sedarweb](https://sedarweb.org/documents/sedar-69-rd13-cookbook-for-using-model-diagnostics-in-integrated-stock-assessments/)
# 
# Large, structured swings in PE are therefore a red flag: they imply that any advice or hindcast relying on the current structural assumptions may have **poor predictive performance and understated risk**.
# 
# ################################################################################
# 
# Ignoring PE diagnostics increases the risk that advice is based on a structurally wrong model with over‑confident uncertainty, which can drive biased catch recommendations and delayed corrective action.
# 
# 
# ## Hidden mis‑specification and biased advice
# 
# - Without checking PE for trends, autocorrelation, or regime shifts, time‑varying productivity, selectivity, or catch reporting changes can be misinterpreted as stable dynamics. [repository.library.noaa](https://repository.library.noaa.gov/view/noaa/52235/noaa_52235_DS1.pdf)
# - This can bias reference points (e.g. \(F_{\text{MSY}}\), \(B_{\text{MSY}}\)) and stock‑status estimates, so catches are set too high when productivity has declined or too low when it has increased, undermining risk equivalence across stocks. [onlinelibrary.wiley](https://onlinelibrary.wiley.com/doi/10.1111/faf.12550)
# 
# ## Over‑confident uncertainty and poor prediction skill
# 
# - Standard likelihoods assume independent, homoscedastic residuals; if PE has low‑frequency structure, variance and autocorrelation are underestimated, so confidence intervals around biomass and forecasts are too narrow. [spo.nmfs.noaa](https://spo.nmfs.noaa.gov/sites/default/files/TMSPO240.pdf)
# - Hindcasting studies show that models with unexamined structured PE often have poor out‑of‑sample prediction skill, leading to optimistic rebuilding trajectories and underestimation of the probability of falling below limit reference points. [sciencedirect](https://www.sciencedirect.com/science/article/abs/pii/S0165783616301540)
# 
# ## Consequences for management outcomes
# 
# - **Higher probability of overfishing and limit breaches**: Biased forecasts combined with narrow uncertainty bands reduce the apparent risk of exceeding \(F_{\text{MSY}}\) or dropping below \(B_{\text{lim}}\), encouraging TACs that are less precautionary than intended. [edepot.wur](https://edepot.wur.nl/556274)
# - **Slow detection of problems**: If retrospective patterns or PE signals are ignored, managers may only react after strong declines are evident in surveys or catches, by which time rebuilding is longer and more costly. [academic.oup](https://academic.oup.com/icesjms/article-pdf/82/2/fsaf014/61741528/fsaf014.pdf)
# - **Inconsistent treatment of data‑limited stocks**: For Category‑2 surplus‑production assessments, failing to use PE diagnostics can make these look more reliable than they are, weakening the precautionary buffer that WKLIFE and ICES seek to maintain for data‑limited advice. [ppl-ai-file-upload.s3.amazonaws](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/11254397/3b640e6f-75f7-479e-aa1f-f7dba5219cd7/wklife_paper_submitted.pdf)
# 
# In short, not using PE diagnostics increases the chance of systematically biased advice, under‑estimated risk, and delayed management response, particularly for data‑limited or structurally complex stocks.
# 
# 
# 
# library(FLRebuild)
# library(r4ss)
# library(FLife)
# 
# source("C:/active/flr/FLRebuild/R/ss-eqPlot.R", echo = TRUE)

compareProdAndPE_plots <- function(tseries,
                                   eqlYield,
                                   shape_pt,
                                   Bmsy_SSB = NULL,
                                   point_col = "grey40",
                                   eq_col = "black",
                                   pt_col_B = "blue",
                                   pt_col_SSB = "red") {
  
  ## =========================
  ## Checks
  ## =========================
  if (!all(c("year", "biomass", "ssb", "yield") %in% names(tseries))) {
    stop("tseries must contain: year, biomass, ssb, yield")
  }
  if (!all(c("ssb", "yield") %in% names(eqlYield))) {
    stop("eqlYield must contain: ssb, yield")
  }
  if (!is.numeric(shape_pt) || length(shape_pt) != 1L ||
      shape_pt <= 0 || shape_pt >= 1) {
    stop("shape_pt must be a single Bmsy/K value in (0,1)")
  }
  if (!is.null(Bmsy_SSB)) {
    if (!is.numeric(Bmsy_SSB) || length(Bmsy_SSB) != 1L || Bmsy_SSB <= 0) {
      stop("Bmsy_SSB must be NULL or a single positive number")
    }
  }
  if (is.null(tseries$pf)) {
    stop("tseries must contain column 'pf' (reference production function values)")
  }
  if (is.null(tseries$pe)) {
    warning("tseries has no 'pe'; reference PE will be recomputed from pf")
  }
  
  ## =========================
  ## Helpers
  ## =========================
  calcSurplusProduction <- function(biomassT, biomassT1, catchT) {
    biomassT1 - biomassT + catchT
  }
  
  ## FLRebuild-style mapping from Bmsy/K to p (allowing p<0 for shape<~0.37)
  getPfromShape <- function(shape) {
    fn <- function(x, y) (y - (1 / (1 + x))^(1 / x))^2
    if (shape < 0.3678794)
      optimise(fn, c(-0.9999, -1e-20), y = shape)$minimum
    else
      optimise(fn, c(1e-20, 10), y = shape)$minimum
  }
  
  ## PT production in r-K-p form (Pella–Tomlinson)
  ptProd_rKp <- function(B, r, K, p) {
    out <- (r / p) * B * (1 - (B / K)^p)
    out[!is.finite(out)] <- NA_real_
    out
  }
  
  ## Biomass PT fit (r, K; p from shape)
  fitPT_biomass <- function(biomassT, surplusProd, shape,
                            weights = NULL, control = list()) {
    
    ok <- is.finite(biomassT) & is.finite(surplusProd) &
      biomassT > 0 & surplusProd > 0
    biomassT    <- biomassT[ok]
    surplusProd <- surplusProd[ok]
    
    if (length(biomassT) < 3L)
      stop("Not enough positive SP data to fit biomass PT")
    
    if (is.null(weights)) {
      weights <- rep(1, length(surplusProd))
    } else {
      weights <- weights[ok]
    }
    
    p <- getPfromShape(shape)
    spLog <- log(surplusProd)
    
    parStart <- c(
      logR = log(max(surplusProd, na.rm = TRUE) /
                   max(biomassT,    na.rm = TRUE)),
      logK = log(max(biomassT,    na.rm = TRUE) * 1.2)
    )
    
    controlDefault <- list(maxit = 1000)
    control <- modifyList(controlDefault, control)
    
    negLogLik <- function(par) {
      r <- exp(par[1])
      K <- exp(par[2])
      
      spHat <- ptProd_rKp(biomassT, r = r, K = K, p = p)
      if (any(!is.finite(spHat))) return(1e12)
      spHat[spHat <= 0] <- .Machine$double.eps
      
      res    <- spLog - log(spHat)
      sigma2 <- sum(weights * res^2) / sum(weights)
      if (!is.finite(sigma2) || sigma2 <= 0) return(1e12)
      
      Kref    <- max(biomassT, na.rm = TRUE)
      penalty <- if (K > 5 * Kref) 5 * (log(K / (5 * Kref)))^2 else 0
      
      0.5 * sum(weights * (log(2 * pi * sigma2) + res^2 / sigma2)) + penalty
    }
    
    fit <- optim(par = parStart,
                 fn  = negLogLik,
                 method  = "Nelder-Mead",
                 control = control,
                 hessian = TRUE)
    
    rHat <- exp(fit$par[1])
    KHat <- exp(fit$par[2])
    p    <- getPfromShape(shape)
    Bmsy <- shape * KHat
    MSY  <- ptProd_rKp(Bmsy, rHat, KHat, p)
    
    list(
      estimates = list(
        r     = rHat,
        k     = KHat,
        p     = p,
        shape = shape,
        bmsy  = Bmsy,
        msy   = MSY
      ),
      convergence = fit$convergence,
      logLik     = -fit$value,
      hessian    = fit$hessian
    )
  }
  
  ## SSB PT fit (r,K or r|K fixed)
  fitPT_ssb <- function(biomassT, surplusProd, shape, Bmsy_input = NULL,
                        weights = NULL, control = list()) {
    
    ok <- is.finite(biomassT) & is.finite(surplusProd) &
      biomassT > 0 & surplusProd > 0
    biomassT    <- biomassT[ok]
    surplusProd <- surplusProd[ok]
    
    if (length(biomassT) < 3L)
      stop("Not enough positive SP data to fit SSB PT")
    
    if (is.null(weights)) {
      weights <- rep(1, length(surplusProd))
    } else {
      weights <- weights[ok]
    }
    
    p     <- getPfromShape(shape)
    spLog <- log(surplusProd)
    
    controlDefault <- list(maxit = 1000)
    control <- modifyList(controlDefault, control)
    
    if (is.null(Bmsy_input)) {
      
      parStart <- c(
        logR = log(max(surplusProd, na.rm = TRUE) /
                     max(biomassT,    na.rm = TRUE)),
        logK = log(max(biomassT,    na.rm = TRUE) * 1.2)
      )
      
      negLogLik <- function(par) {
        r <- exp(par[1])
        K <- exp(par[2])
        
        spHat <- ptProd_rKp(biomassT, r = r, K = K, p = p)
        if (any(!is.finite(spHat))) return(1e12)
        spHat[spHat <= 0] <- .Machine$double.eps
        
        res    <- spLog - log(spHat)
        sigma2 <- sum(weights * res^2) / sum(weights)
        if (!is.finite(sigma2) || sigma2 <= 0) return(1e12)
        
        Kref    <- max(biomassT, na.rm = TRUE)
        penalty <- if (K > 5 * Kref) 5 * (log(K / (5 * Kref)))^2 else 0
        
        0.5 * sum(weights * (log(2 * pi * sigma2) + res^2 / sigma2)) + penalty
      }
      
      fit <- optim(par = parStart,
                   fn  = negLogLik,
                   method  = "Nelder-Mead",
                   control = control,
                   hessian = TRUE)
      
      rHat <- exp(fit$par[1])
      KHat <- exp(fit$par[2])
      
    } else {
      
      KHat <- Bmsy_input / shape
      
      parStart <- c(
        logR = log(max(surplusProd, na.rm = TRUE) /
                     max(biomassT,    na.rm = TRUE))
      )
      
      negLogLik <- function(par) {
        r <- exp(par[1])
        
        spHat <- ptProd_rKp(biomassT, r = r, K = KHat, p = p)
        if (any(!is.finite(spHat))) return(1e12)
        spHat[spHat <= 0] <- .Machine$double.eps
        
        res    <- spLog - log(spHat)
        sigma2 <- sum(weights * res^2) / sum(weights)
        if (!is.finite(sigma2) || sigma2 <= 0) return(1e12)
        
        0.5 * sum(weights * (log(2 * pi * sigma2) + res^2 / sigma2))
      }
      
      fit <- optim(par = parStart,
                   fn  = negLogLik,
                   method  = "Nelder-Mead",
                   control = control,
                   hessian = TRUE)
      
      rHat <- exp(fit$par[1])
    }
    
    Bmsy <- shape * KHat
    MSY  <- ptProd_rKp(Bmsy, rHat, KHat, p)
    
    list(
      estimates = list(
        r     = rHat,
        k     = KHat,
        p     = p,
        shape = shape,
        bmsy  = Bmsy,
        msy   = MSY
      ),
      convergence = fit$convergence,
      logLik     = -fit$value,
      hessian    = fit$hessian
    )
  }
  
  predictPtBiomass <- function(biomassT, catchT, r, K, p) {
    biomassT + ptProd_rKp(biomassT, r, K, p) - catchT
  }
  
  ## =========================
  ## Data prep
  ## =========================
  year <- tseries$year
  B    <- tseries$biomass
  SSB  <- tseries$ssb
  C    <- tseries$yield
  
  eq_SSB <- eqlYield$ssb
  eq_SP  <- eqlYield$yield
  
  yr_t   <- year[-length(year)]
  B_t    <- B[-length(B)]
  B_t1   <- B[-1]
  SSB_t  <- SSB[-length(SSB)]
  SSB_t1 <- SSB[-1]
  C_t    <- C[-length(C)]
  
  SP_B   <- calcSurplusProduction(B_t,   B_t1,   C_t)
  SP_SSB <- calcSurplusProduction(SSB_t, SSB_t1, C_t)
  
  ## =========================
  ## PT fits
  ## =========================
  fit_B   <- fitPT_biomass(B_t, SP_B, shape = shape_pt)
  fit_SSB <- fitPT_ssb  (SSB_t, SP_SSB, shape = shape_pt, Bmsy_input = Bmsy_SSB)
  
  est_B   <- fit_B$estimates
  est_SSB <- fit_SSB$estimates
  
  ## =========================
  ## Curves
  ## =========================
  B_grid   <- seq(0, est_B$k,   length.out = 400)
  SSB_grid <- seq(0, est_SSB$k, length.out = 400)
  
  PT_B_grid   <- ptProd_rKp(B_grid,   est_B$r,   est_B$k,   est_B$p)
  PT_SSB_grid <- ptProd_rKp(SSB_grid, est_SSB$r, est_SSB$k, est_SSB$p)
  
  PT_B_grid_plot   <- pmax(PT_B_grid,   0)
  PT_SSB_grid_plot <- pmax(PT_SSB_grid, 0)
  
  eq_keep <- is.finite(eq_SSB) & is.finite(eq_SP) & eq_SP >= 0
  
  ylim_B   <- c(0, max(c(SP_B[SP_B >= 0], PT_B_grid_plot),   na.rm = TRUE) * 1.05)
  ylim_SSB <- c(0, max(c(SP_SSB[SP_SSB >= 0],
                         eq_SP[eq_keep],
                         PT_SSB_grid_plot), na.rm = TRUE) * 1.05)
  
  ## =========================
  ## Predictions & PE
  ## =========================
  Bpred_PT   <- predictPtBiomass(B_t,   C_t, est_B$r,   est_B$k,   est_B$p)
  SSBpred_PT <- predictPtBiomass(SSB_t, C_t, est_SSB$r, est_SSB$k, est_SSB$p)
  
  eps_B_PT <- log(pmax(B_t1,   .Machine$double.eps)) -
    log(pmax(Bpred_PT,   .Machine$double.eps))
  
  eps_SSB_PT <- log(pmax(SSB_t1, .Machine$double.eps)) -
    log(pmax(SSBpred_PT, .Machine$double.eps))
  
  ## Reference PF-based PE (mis-specification diagnostic)
  pf_t   <- tseries$pf[-length(tseries$pf)]
  SSB_t0 <- tseries$ssb[-length(tseries$ssb)]
  SSB_t1_0 <- tseries$ssb[-1]
  C_t0  <- tseries$yield[-length(tseries$yield)]
  
  pe_from_pf <- log(pmax(SSB_t1_0, .Machine$double.eps)) -
    log(pmax(SSB_t0 - C_t0 + pf_t, .Machine$double.eps))
  
  ## If tseries$pe exists, keep it as given
  pe_ref <- if (!is.null(tseries$pe)) tail(tseries$pe, -1) else pe_from_pf
  
  ylim_PE <- range(c(eps_B_PT, eps_SSB_PT, pe_ref), na.rm = TRUE)
  if (!all(is.finite(ylim_PE))) ylim_PE <- c(-1, 1)
  
  ## =========================
  ## Plots
  ## =========================
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))
  par(mar = c(4, 4, 2, 1))
  
  ## Biomass SP
  plot(B_t, SP_B,
       pch = 16, col = point_col,
       xlab = "Biomass", ylab = "Surplus production",
       xlim = c(0, est_B$k), ylim = ylim_B,
       main = paste0("Biomass SP; Bmsy/K=", round(est_B$shape, 2)))
  lines(B_t, SP_B, lwd = 0.25)
  lines(B_grid, PT_B_grid_plot, lwd = 2, col = pt_col_B)
  abline(h = 0, lty = 3)
  
  ## SSB SP
  plot(SSB_t, SP_SSB,
       pch = 16, col = point_col,
       xlab = "SSB", ylab = "Surplus production",
       xlim = c(0, est_SSB$k), ylim = ylim_SSB,
       main = paste0("SSB SP; Bmsy/K=", round(est_SSB$shape, 2)))
  lines(SSB_t, SP_SSB, lwd = 0.25)
  lines(eq_SSB[eq_keep], eq_SP[eq_keep], lwd = 2, col = "orange")
  lines(SSB_grid, PT_SSB_grid_plot, lwd = 2, col = pt_col_SSB)
  abline(h = 0, lty = 3)
  
  ## Process error panel
  plot(yr_t, pe_ref,
       type = "l", lwd = 2, col = "orange",
       xlab = "Year", ylab = "Log process error",
       ylim = ylim_PE,
       main = "Process error (reference vs PT)")
  lines(yr_t, eps_B_PT,   lwd = 2, col = pt_col_B)
  lines(yr_t, eps_SSB_PT, lwd = 2, col = pt_col_SSB)
  abline(h = 0, lty = 2)
  legend("topleft",
         legend = c("Reference PE (pf)", "Biomass PT PE", "SSB PT PE"),
         col    = c("orange", pt_col_B,   pt_col_SSB),
         lwd    = 2,
         lty    = c(1, 1, 1),
         bty    = "n")
  
  ## =========================
  ## Return object
  ## =========================
  pars_df <- data.frame(
    space = c("biomass", "ssb"),
    r     = c(est_B$r,    est_SSB$r),
    K     = c(est_B$k,    est_SSB$k),
    p     = c(est_B$p,    est_SSB$p),
    Bmsy  = c(est_B$bmsy, est_SSB$bmsy),
    MSY   = c(est_B$msy,  est_SSB$msy),
    BmsyK = c(est_B$shape, est_SSB$shape),
    stringsAsFactors = FALSE
  )
  
  pe_df <- data.frame(
    year        = yr_t,
    pe_ref      = pe_ref,
    pe_from_pf  = pe_from_pf,
    pe_B_PT     = eps_B_PT,
    pe_SSB_PT   = eps_SSB_PT
  )
  
  invisible(list(
    pars    = pars_df,
    pe      = pe_df,
    biomass = list(
      fit      = fit_B,
      observed = data.frame(year = yr_t,
                            biomass = B_t,
                            sp      = SP_B,
                            pe_pt   = eps_B_PT),
      curve    = data.frame(biomass = B_grid,
                            pt      = PT_B_grid)
    ),
    ssb = list(
      fit      = fit_SSB,
      observed = data.frame(year = yr_t,
                            ssb     = SSB_t,
                            sp      = SP_SSB,
                            pe_pt   = eps_SSB_PT),
      curve    = data.frame(ssb = SSB_grid,
                            pt  = PT_SSB_grid)
    )
  ))
}

if (FALSE){
  baSS=SS_output("C:/active/flrpapers/PE/SS3/bass-north")
  bass=curveSS(baSS)
  bass$curve=bass$curve[order(bass$curve$ssb),]
  
  bass$curve$ssb[which.max(bass$curve$yield)]/ max(bass$curve$ssb)
  
  res=compareProdAndPE_plots(bass$tseries,bass$curve,Bmsy_SSB=bass[[3]]$bmsy,shape_pt=0.2)
  }  
