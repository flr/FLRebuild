#' Build equilibrium and production-curve inputs from SS output
#'
#' Process Stock Synthesis output to calculate yield and surplus production
#' series for equilibrium-curve plotting.
#'
#' @param object A list returned by `r4ss::SS_output()`, or a list of character
#'   paths to Stock Synthesis run directories.
#' @param maxY Multiplier used to set the upper y-limit proxy for triangle data.
#'
#' @return A list with elements `tseries`, `curve`, `refpts`, `triangle`, and
#'   `derived`. For list input containing paths, each data.frame includes an
#'   identifying column `id`.
#' @export
setMethod("curveSS", signature(object = "list"), function(object, maxY = 1.5) {
  # If input is a list of character paths, dispatch per element and bind with ID.
  if (length(object) > 0 && all(vapply(object, function(x) is.character(x) && length(x) == 1, logical(1)))) {
    ids <- names(object)
    if (is.null(ids) || any(ids == "")) ids <- as.character(seq_along(object))

    runs <- lapply(seq_along(object), function(i) {
      out <- curveSS(object[[i]], maxY = maxY)
      lapply(out, function(df) {
        df <- as.data.frame(df)
        cbind(id = ids[i], df, stringsAsFactors = FALSE)
      })
    })

    parts <- names(runs[[1]])
    combined <- lapply(parts, function(p) do.call(rbind, lapply(runs, `[[`, p)))
    names(combined) <- parts
    return(combined)
  }

  eqlYield  =object[["equil_yield"]]
  timeseries=object[["timeseries"]]
  rfs=c("SSB_unfished","Totbio_unfished","SmryBio_unfished","Recr_unfished","SSB_Btgt","SPR_Btgt","annF_Btgt",  
        "Dead_Catch_Btgt","SSB_SPR", "annF_SPR","Dead_Catch_SPR","SSB_MSY","SPR_MSY","annF_MSY",   
        "Dead_Catch_MSY", "Ret_Catch_MSY","B_MSY/SSB_unfished")
  dq=object[["derived_quants"]][unique(object$derived_quants$Label)[rfs],]
  
  if (!("Tot_Catch"%in%names(eqlYield))&("Catch"%in%names(eqlYield)))
    names(eqlYield)[names(eqlYield)=="Catch"]="Tot_Catch"

  if (any(tolower(names(timeseries))%in%"yr"))
    names(timeseries)[tolower(names(timeseries))=="yr"]="year"

  # Filter and prepare data for plotting
  ts <- timeseries %>%
    dplyr::filter(!Era %in% c("VIRG", "FORE")) %>%
    dplyr::mutate(catch_total = rowSums(dplyr::select(., dplyr::starts_with("dead(B)")), na.rm = TRUE)) %>%
    dplyr::group_by(year, Seas) %>%
    dplyr::summarise(
      sum_bio_all = sum(Bio_all, na.rm = TRUE),
      sum_spawn_bio = sum(SpawnBio, na.rm = TRUE),
      sum_catch_total = sum(catch_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      mean_bio_all = mean(sum_bio_all),
      mean_spawn_bio = mean(sum_spawn_bio),
      catch_total = sum(sum_catch_total),
      .groups = "drop"
    )
  names(ts)=c("year","biomass","ssb","yield")
  
  nyears=nrow(ts)
  
  ts$sp   =c(NA, ts$biomass[-1] - ts$biomass[-nyears] + ts$yield[-nyears])
  eql     =dplyr::transmute(eqlYield, yield=Tot_Catch,ssb=SSB)
  rfs     =dplyr::transmute(subset(eqlYield, Tot_Catch==max(Tot_Catch)), msy=Tot_Catch,bmsy=SSB)[1,]
  maxY    =signif(max(c(ts$yield,eql$yield))*maxY,1)
  triangle=data.frame(x=c(rfs$bmsy, rfs$bmsy, rfs$bmsy*maxY/rfs$msy, rfs$bmsy),
                      y=c(rfs$msy,      maxY,                  maxY, rfs$msy))
  vBiomass=vBio(object)
  
  prodFun<-splinefun(x=eql$ssb, y =eql$yield, method = "natural")
  
  ts=cbind(ts,pf=prodFun(ts$ssb))
  ts=cbind(ts,pe=c(with(ts, (ssb[-1]-tail(ssb,-1)-tail(yield,-1)+tail(pf,-1))/tail(ssb,-1)),NA))

  return(list(tseries=merge(ts,vBiomass),curve=eql,refpts=rfs,triangle=triangle,derived=dq))
})

#' Build equilibrium and production-curve inputs from SS output directory
#'
#' Read Stock Synthesis output from a directory and return the same structure as
#' `curveSS` for a list object.
#'
#' @param object Character path (or vector of paths) to Stock Synthesis run
#'   directories.
#' @param maxY Multiplier used to set the upper y-limit proxy for triangle data.
#' @param covar Logical; passed to `r4ss::SS_output`.
#'
#' @return A list with elements `tseries`, `curve`, `refpts`, `triangle`, and
#'   `derived`.
#' @export
setMethod("curveSS", signature(object = "character"), function(object, maxY = 1.5, covar = FALSE) {
  if (length(object) > 1) {
    ids <- names(object)
    if (is.null(ids) || any(ids == "")) ids <- as.character(seq_along(object))

    runs <- lapply(seq_along(object), function(i) {
      out <- curveSS(object[[i]], maxY = maxY, covar = covar)
      lapply(out, function(df) {
        df <- as.data.frame(df)
        cbind(id = ids[i], df, stringsAsFactors = FALSE)
      })
    })

    parts <- names(runs[[1]])
    combined <- lapply(parts, function(p) do.call(rbind, lapply(runs, `[[`, p)))
    names(combined) <- parts
    return(combined)
  }

  ss_rep <- r4ss::SS_output(dir = object, covar = covar, verbose = FALSE, printstats = FALSE)
  curveSS(ss_rep, maxY = maxY)
})



logistic <- plogis  # maps to standard logistic function
logit <- qlogis     # maps to inverse logistic function

#' Estimate vulnerable biomass from SS report components
#'
#' Computes sex-specific vulnerable biomass using selectivity parameters and
#' summary biomass time series from Stock Synthesis output.
#'
#' @param object A list returned by `r4ss::SS_output()`.
#'
#' @return A data frame with columns `year`, `female`, and `male`.
#' @export
setMethod("vBio", signature(object = "list"), function(object) {
  # Get selectivity parameters
  sel_F_peak  =unlist(c(subset(object$parameters, Label=="SzSel_Fem_Peak_Fishery_2(2)","Value")))
  sel_F_ascend=unlist(c(subset(object$parameters, Label=="SzSel_Fem_Ascend_Fishery_2(2)","Value")))
  sel_M_peak  =unlist(c(subset(object$parameters, Label=="SzSel_Male_Peak_Fishery_3(3)","Value")))
  sel_M_ascend=unlist(c(subset(object$parameters, Label=="SzSel_Male_Ascend_Fishery_3(3)","Value")))
  
  # Calculate selectivity using logistic transform
  calc_sel<-function(len, peak, ascend) 
    logistic((len - peak)/exp(ascend))
  
  # Get length bins from data
  len_bins=as.numeric(gsub("\\D", "", names(object$sizeselex)[-(1:4)]))
  len_bins[is.na(len_bins)]=0
  
  # Calculate sex-specific selectivity
  sel_F=calc_sel(len_bins, sel_F_peak, sel_F_ascend)
  sel_M=calc_sel(len_bins, sel_M_peak, sel_M_ascend)
  
  # Extract and scale biomass (assuming SmryBio is in log space)
  bio_F=object$timeseries$`SmryBio_SX:1_GP:1`
  bio_M=object$timeseries$`SmryBio_SX:2_GP:1`
  
  # Calculate vulnerable biomass
  rtn=data.frame(
    year   = object$timeseries$Yr,
    female = bio_F * sum(sel_F),
    male   =NA)
  
  if(!is.null(bio_M)&!is.null(sel_M))
    rtn$male = bio_M * sum(sel_M)

  return(rtn)
})
