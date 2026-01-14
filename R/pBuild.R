#' Predict Rebuilding Times from SSB Levels
#'
#' @description
#' Interpolates rebuilding times (time to reach target SSB) for given initial SSB levels
#' by using rebuild trajectory results. The function calculates recovery times from rebuild
#' trajectories and then interpolates these times for the provided SSB values.
#'
#' @param ssb An FLQuant object containing spawning stock biomass values for which
#'   rebuilding times are to be predicted
#' @param eq An FLBRP object containing the equilibrium/biological reference points
#'   used to generate rebuild trajectories
#' @param ftar Numeric scalar; target F as a multiple of MSY harvest used when
#'   computing rebuild trajectories. Default is 0 (F = 0).
#' @param rb Optional precomputed rebuild object. By default
#'   \code{rebuild(eq, targetF = ftar * refpts(eq)["msy", "harvest"])} is used.
#'
#' @return An FLQuant object with the same dimensions as the input `ssb`, containing
#'   the predicted rebuilding times (time to reach target SSB, relative to MSY SSB).
#'   Values are set to 0 if rebuilding time cannot be determined.
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Calling \code{rebuild(eq)} to generate rebuilding trajectories from the equilibrium object
#'   \item Calculating initial SSB levels (relative to MSY SSB) and corresponding recovery times
#'   \item Using linear interpolation to predict rebuilding times for the provided SSB values
#'   \item Returning an FLQuant with predicted recovery times matching the input SSB dimensions
#' }
#'
#' The function handles errors gracefully, returning an FLQuant of zeros with matching
#' dimensions if interpolation cannot be performed (e.g., insufficient data or invalid inputs).
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' library(FLBRP)
#'
#' # Create or load an FLBRP object
#' # eq <- brp(flstock_object)
#'
#' # Create SSB values for prediction
#' # ssb_values <- FLQuant(seq(0.1, 0.9, by=0.1) * refpts(eq)["msy","ssb"])
#'
#' # Predict rebuilding times
#' # rebuild_times <- pBuild(ssb_values, eq)
#' }
#'
#' @seealso \code{\link{rebuild}} for generating rebuilding trajectories
#' @importFrom plyr ddply
#' @importFrom data.table as.data.table
#' @importFrom FLCore as.FLQuant
#' @export
pBuild <- function(ssb, eq, 
                   ftar =0, 
                   rb   =rebuild(eq, targetF = ftar * refpts(eq)["msy", "harvest"]),
                   btar =refpts(eq)["msy","ssb"]) {

  ## Recovery to (data.y) by year from initial biomass (data.x)
  t1=as.data.frame(ssb(rb)[,1]%/%btar)
  t2=as.data.frame(ssb(rb)    %/%btar)
  t3=merge(t1[,-2],t2,by=names(t1)[c(1,3:6)])

  ## Relationship between initial biomass and recovery level by year
  t4=plyr::ddply(t3,.(data.x), with, 
                  suppressWarnings(tryIt(approx(data.y,year,1)$y-1)))
  names(t4)=c("initial","tRecover")
    
  # Relative state
  t5=ssb%/%btar
  t5=as.data.frame(t5)
  
  t4=t4[complete.cases(t4),]
  t5=t5[complete.cases(t5),]
  t5=t5[t5$data>min(t4$initial)|t5$data<max(t4$initial),]
  
  estimateRecovery<-function(t4, t5) {
    # Working with data.table
    t5_dt=data.table::as.data.table(t5)
    
    # Extract the interpolation grid
    x0=t4$initial
    y0=t4$tRecover
    
    # for each (iter, year) group, interpolate over that group's "data"
    t5_dt[, .(data = approx(x=x0, y=y0, xout=data)$y), by=.(iter, year)]}
  
  rtn=tryIt(estimateRecovery(t4, t5))
  
  if (is.null(rtn)) 
    return(FLQuant()) 
  
  # Convert data.table to data.frame then to FLQuant
  rtnDF=data.frame(year=rtn$year,
                    iter=rtn$iter,
                    data=rtn$data)
  
  rtn = FLCore::as.FLQuant(rtnDF)

  if (any(ssb%/%btar>1))
      rtn[ssb%/%btar>1]=0
  
  rtn}
