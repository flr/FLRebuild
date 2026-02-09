#' Calculate Leslie Matrix Demographic Properties
#'
#' This function computes demographic parameters from a Leslie population matrix,
#' including damping ratio, intrinsic growth rate, generation time, and net reproductive rate.
#'
#' @param x An FLBRP object containing life history information
#' @param fbar An FLQuant object specifying the fishing mortality rate. 
#'   Default is MSY harvest rate from reference points.
#'
#' @return A named numeric vector containing:
#'   \item{damp}{Damping ratio - the ratio of the first to second largest eigenvalue moduli.
#'     Indicates how quickly a population returns to stable age distribution after perturbation.}
#'   \item{r}{Intrinsic population growth rate (log of dominant eigenvalue)}
#'   \item{gt}{Generation time - mean age of reproduction weighted by survivorship and fecundity}
#'   \item{r0}{Net reproductive rate (R0) - expected lifetime reproductive output}
#'
#' @details
#' The function constructs a Leslie matrix using FLife::leslie() and derives key demographic
#' parameters. The damping ratio measures resilience to age structure perturbations. Generation
#' time represents the average time between generations weighted by reproductive contributions.
#'
#' @examples
#' \dontrun{
#' # Calculate demographic properties at MSY
#' demo_params <- leslieFn(brp_object)
#' 
#' # Calculate at a specific fishing mortality
#' demo_params <- leslieFn(brp_object, fbar = FLQuant(0.2))
#' }
#'
#' @export
leslieFn<-function(x,fbar=FLQuant(c(refpts(x)["msy","harvest"]))){
  
  L=FLife::leslie(x,fbar=fbar)[drop=TRUE]
  
  eigvals=eigen(L)$values
  
  egnMod=Mod(eigvals)
  
  sorted=sort(egnMod, decreasing=TRUE)
  if(length(sorted)>1)
    dampingRatio=sorted[1]/sorted[2]
  else 
    dampingRatio=NA
  
  r=log(lambda(L))
  
  
  fecundity=L[1, ]
  ages     =length(fecundity)
  

  lx   =numeric(ages)
  lx[1]=1
  for(i in 2:ages)
    lx[i]=lx[i-1] * L[i, i-1]
  
  # Net reproductive rate (R0)
  R0=sum(lx * fecundity)
  
  # Generation time (mean age of reproduction)
  genTime=sum((1:ages) * lx * fecundity) / R0
  
  c(damp=dampingRatio,r=r,gt=genTime,r0=R0)}


#L=FLife::leslie(bh1[[1]],fbar=FLQuant(c(refpts(bh1[[1]])["msy","harvest"]))[drop=TRUE])

#' Calculate Life History and Demographic Covariates
#'
#' This function computes a comprehensive set of life history traits and demographic
#' parameters from an FLBRP object, useful for comparative analyses and building
#' statistical models of population dynamics and rebuilding capacity.
#'
#' @param object An FLBRP object containing stock assessment information and life history data
#' @param fbar An FLQuant specifying fishing mortality. Default is MSY harvest rate from reference points.
#' @param model Character string specifying the stock-recruitment model type. 
#'   Options: "bevholt" (default) or "ricker". Used for calculating steepness.
#'
#' @return A named numeric vector containing life history and demographic covariates:
#'   
#' **Weight and Age Metrics:**
#'   \item{swt}{Mean spawner weight-at-age (weighted by spawning stock numbers and maturity)}
#'   \item{cwt}{Mean catch weight (weighted by catch numbers)}
#'   \item{ssbAge}{Age at which spawning stock biomass contribution is maximized (apex age)}
#'   \item{ctcAge}{Age at which catch biomass contribution is maximized}
#'   
#' **Timing and Selectivity:**
#'   \item{m.spwn}{Logical - whether natural mortality before spawning is non-zero}
#'   \item{harvest.spwn}{Logical - whether fishing mortality before spawning is non-zero}
#'   \item{wts}{Mean absolute relative difference between catch and stock weights}
#'   \item{ks}{Kolmogorov-Smirnov test statistic comparing selectivity and maturity schedules}
#'   \item{matsel}{Area between selectivity and maturity curves (using ECDF comparison)}
#'   \item{abiSSB}{Area between ABI (aggregate biomass index) and SSB age distributions}
#'   
#' **Productivity and Reference Points:**
#'   \item{spr0}{Spawning potential ratio at unfished conditions}
#'   \item{shape}{Production function shape (ratio of B_MSY to virgin biomass)}
#'   \item{s}{Steepness parameter of the stock-recruitment relationship}
#'   \item{pe}{Process error - standard deviation of recruitment residuals. Measures variability in stock-recruitment relationship}
#'   
#' **Demographic Parameters (from leslieFn):**
#'   \item{damp}{Damping ratio from Leslie matrix eigenvalues}
#'   \item{r}{Intrinsic population growth rate}
#'   \item{gt}{Generation time (mean age of reproduction)}
#'   \item{r0}{Net reproductive rate}
#'
#' @details
#' The function performs multiple calculations:
#' 
#' 1. **Apex ages**: Identifies ages contributing most to catch and SSB
#' 2. **Mean weights**: Calculates abundance-weighted mean weights for spawners and catch
#' 3. **Selectivity-maturity comparison**: Uses ECDF (empirical cumulative distribution function)
#'    to quantify differences between selectivity and maturity patterns
#' 4. **ABI-SSB comparison**: Compares age distributions based on aggregate biomass index
#'    versus spawning stock biomass contributions
#' 5. **Demographic properties**: Calls leslieFn() to compute Leslie matrix-based metrics
#'
#' The metrics are useful for:
#' - Comparative life history analysis across species/stocks
#' - Building statistical models of rebuild time and recovery capacity
#' - Understanding relationships between life history and productivity
#' - Identifying key biological traits affecting management outcomes
#' - Quantifying uncertainty in stock-recruitment relationships via process error
#'
#' @examples
#' \dontrun{
#' # Calculate covariates at MSY fishing mortality
#' covars <- covarFn(brp_object)
#' 
#' # Calculate for Ricker stock-recruitment relationship
#' covars <- covarFn(brp_object, model = "ricker")
#' 
#' # Calculate at different fishing mortality
#' covars <- covarFn(brp_object, fbar = as.FLQuant(0.3))
#' 
#' # Use in comparative analysis across multiple stocks
#' covar_list <- lapply(stock_list, covarFn)
#' covar_df <- do.call(rbind, covar_list)
#' }
#'
#' @seealso 
#' \code{\link{leslieFn}} for Leslie matrix demographic calculations
#' 
#' @references
#' Kolmogorov-Smirnov test for comparing distributions
#' 
#' @export
covarFn<-function(object,fbar=as.FLQuant(refpts(object)["msy","harvest",drop=T]),model="bevholt"){
  
  smry<-function(x){
    
    # Mat & Sel
    ecdf1    =ecdf(c(catch.sel(x)/max(catch.sel(x))))
    ecdf2    =ecdf(c(mat(x)))
    rangeVals=seq(min(c(catch.sel(x)/max(catch.sel(x)),mat(x))),
                  max(c(catch.sel(x)/max(catch.sel(x)),mat(x))),
                  length.out=1000)
    
    # Evaluate ECDFs over the range
    ecdfDiff=ecdf1(rangeVals)-ecdf2(rangeVals)
    # Area Between the Curves
    matsel=sum(ecdfDiff)*diff(rangeVals)[1]
    
    # ABI & SSB
    abiSel=as.numeric(c(ages(catch.sel(x)))>=c(rebuild::abiAge(x)))
    ssbSel=c(mat(x)%*%stock.wt(x))
    ssbSel=ssbSel/max(ssbSel)
    
    ecdf1 =ecdf(abiSel)
    ecdf2 =ecdf(ssbSel)
    rangeVals=seq(0,1,length.out=100)
    
    # Evaluate ECDFs over the range
    ecdfDiff=ecdf1(rangeVals)-ecdf2(rangeVals)
    # Area Between the Curves
    abiSSB=sum(ecdfDiff)*diff(rangeVals)[1]
    
    rtn=c(m.spwn      =any(m.spwn(x)!=0),
          harvest.spwn=any(harvest.spwn(x)!=0),
          wts         =mean(catch.sel(x)%*%(catch.wt(x)-stock.wt(x))/stock.wt(x)),
          ks          =ks.test(catch.sel(x)/max(catch.sel(x)),mat(x))[[1]],
          matsel      =matsel,
          abiSSB      =abiSSB,
          spr0        =spr0(x),
          shape       =refpts(x)["msy","ssb"]/refpts(x)["virgin","ssb"],
          s           =tryIt(sv(params(x),spr0=spr0(x),model=model)["s"]),
          pe          =var(rebuild::processError(x)$pe,na.rm=TRUE)^0.5)
        
    rtn=c(rtn,rebuild::leslieFn(x))
    names(rtn)=gsub("leslieFn.","",names(rtn))
    
    rtn}
  
  ssbN<-function(object) {
    
    f = FLCore::harvest(object) %*% harvest.spwn(object)
    m = m(object) %*% m.spwn(object)
    
    expZ = exp(-f %-% m)
    
    res = (stock.n(object) %*% expZ %*% stock.wt(object) %*%
             mat(object))
    
    return(res)}
  
  fbar(object)=fbar
  
  ags=ages(stock.n(object))
  
  apexCtc  =catch.n(object)%*%catch.wt(object)
  apexCtc[]=fapex(apexCtc)
  apexCtc=apexCtc==catch.n(object)%*%catch.wt(object)
  apexCtc=subset(as.data.frame(as.FLQuant(apexCtc)),data==1)[,-7]
  
  apexSSB  =ssbN(object)
  apexSSB[]=fapex(ssbN(object))
  apexSSB  =apexSSB==ssbN(object)
  apexSSB=subset(as.data.frame(as.FLQuant(apexSSB)),data==1)[,-7]
  
  names(apexCtc)[1]="ctcAge"
  names(apexSSB)[1]="ssbAge"
  
  mnCwt=as.data.frame(quantSums(catch.n(object)%*%catch.wt(object))%/%
                        quantSums(catch.n(object)))[,-1]
  mnSwt=as.data.frame(quantSums(stock.n(object)%*%stock.wt(object)%*%mat(object))%/%
                        quantSums(stock.n(object)%*%mat(object)))[,-1]
  
  names(mnSwt)[6]="swt"
  names(mnCwt)[6]="cwt"
  
  rtn=unlist(merge(merge(mnSwt,  mnCwt,  by=names(mnSwt)[  1:5]),
                   merge(apexSSB,apexCtc,by=names(apexSSB)[2:6]),by=names(mnSwt)[1:5])[,6:9])
  
  
  c(rtn,smry(object))}

if(FALSE){
  
  load("C:/active/flrpapers/rebuild/papers/brebuild/data/age/eql.RData")
  
  covars=rbind(
    cbind(SRR="Beverton & Holt",rbind(
      cbind(Scenario=1,ldply(bh1, function(x) tryIt(covarFn(x)))),
      cbind(Scenario=2,ldply(bh2, function(x) tryIt(covarFn(x)))),
      cbind(Scenario=3,ldply(bh3, function(x) tryIt(covarFn(x)))))),
    cbind(SRR="Ricker",rbind(
      cbind(Scenario=1,ldply(rk1, function(x) tryIt(covarFn(x,model="ricker")))),
      cbind(Scenario=2,ldply(rk2, function(x) tryIt(covarFn(x,model="ricker")))),
      cbind(Scenario=3,ldply(rk3, function(x) tryIt(covarFn(x,model="ricker")))))))
  covars$Scenarios=factor(covars$Scenario, 
                          labels=c("h & R0 Estimated", "Steepness Prior", "R0 Fixed"))
  
  save(covars,file="C:/active/flrpapers/rebuild/papers/brebuild/data/age/covars.RData") 
}