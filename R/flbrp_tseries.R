#' @title tseries
#' 
#' @description Calculates the surplus production and expected yield etc for the estimates of SSB and biomass
#'
#' @param object an \code{FLBRP} object 
#' @param seasons a numeric with seasons
#' 
#' @aliases
#' 
#' @return \code{FLQuants} object
#'
#' @seealso \code{\link{expand}}
#'
#' @importFrom plyr laply alply ldply
#' @importFrom FLCore mcf
#' @export tseries
#' @docType methods
#' @rdname tseries
#'
#' 
#' @examples
#' \dontrun{
#' }

setMethod("tseries", signature(object="FLBRP"), function(object){
  ebiomass.obs<-function(x) attributes(x)$eb.obs

  nms=dimnames(refpts(object))
  nms$refpt=paste("ssb",dimnames(ssb.obs(object))$year,sep="")
  
  vrgn=refpts(object)["virgin","ssb"]
  
  discards.obs(object)[is.na(discards.obs(object))]=0
  
  landings.sel(object)=landings.sel(object)+discards.sel(object)
  discards.sel(object)[]=0
  
  rfs=FLPar(array(NA,laply(nms,length),dimnames=nms))
  rfs[,"ssb",]=ssb.obs(object)
  refpts(object)=rfs
  rtn=refptsEB(object)
  fbar(object)=FLQuant(c(rtn[,"harvest"]))
  
  rtn=alply(rtn,2,FLQuant,dimnames=dimnames(ssb.obs(object)))
  names(rtn)=as.character(unlist(attributes(rtn)$split_labels))
  
  rtn$eb =ebiomass.obs(object)
  
  rtn$spSSB=ssb.obs(object)[,-1]-ssb.obs(object)[,-dim(ssb.obs(object))[2]]+catch.obs(object)[,-dim(ssb.obs(object))[2]]
  rtn$spEB =ebiomass.obs(object)[,-1]-ebiomass.obs(object)[,-dim(ebiomass.obs(object))[2]]+catch.obs(object)[,-dim(ebiomass.obs(object))[2]]
  
  rtn=mcf(as(rtn,"FLQuants"))
  
  ord=c("harvest","yield","rec","ssb","eb","revenue","cost","profit","spSSB","spEB")
  
  rtn=rtn[ord]
  rtn[["peSSB"]]=rtn$spSSB-rtn$yield
  rtn[["peSSB"]]=rtn[["peSSB"]]%/%rtn[["ssb"]]
  rtn[["peEB"]] =rtn$spEB -rtn$yield
  rtn[["peEB"]]=rtn[["peEB"]]%/%rtn[["eb"]]
  
  chk=ssb.obs(object)
  chk=chk%=%vrgn
  chk=chk<=ssb.obs(object)
  
  rtn[["ssb"]]=ssb.obs(object)
  
  if(any(chk)){
    for (i in names(rtn))
      rtn[[i]][chk]=NA}
  
  rtn})

setMethod("tseries", signature(object="FLBRPs"), function(object){
      ldply(object, function(x) model.frame(tseries(x)))})
   
setGeneric("prodPts", function(object, ...) standardGeneric("prodPts"))
setMethod( "prodPts", signature(object="FLBRP"),  function(object){ tseries(object)})
setMethod( "prodPts", signature(object="FLBRPs"), function(object){
  ldply(object, function(x) model.frame(prodPts(x)))})

setGeneric("prodFn", function(object, ...) standardGeneric("prodFn"))
setMethod( "prodFn", signature(object="FLBRPs"), function(object){
  ldply(object, function(x) model.frame(prodFn(x)))})

setMethod("prodFn",signature(object="FLBRP"),function(object) {
  
  # Extract data
  pf=model.frame(FLQuants(object,
                          "f"    =fbar, 
                          "yield"=function(x) catch(x)))
  maxF=subset(pf,f>0&yield==0)[1,"f"]
  if (!is.na((maxF)))
    fbar(object)=FLQuant(seq(0,maxF,length.out=101))

  pf=model.frame(FLQuants(object,
                          "harvest" =fbar, 
                          "yield"   =function(x) catch(x),
                          "rec"     =function(x) rec(x),
                          "ssb"     =function(x) ssb(x),
                          "biomass" =function(x) stock(x),
                          "eb"      =function(x) ebiomass(x)))
  
  return(pf)})

setGeneric("MLP", function(object,...) standardGeneric("MLP"))
setMethod("MLP", signature(object="FLBRPs"), function(object,rm.neg=TRUE){
  ldply(object, function(x) model.frame(MLP(x,rm.neg)))})

setMethod("MLP",signature(object="FLBRP"), 
    function(object,rm.neg=TRUE) {
            
        pt=model.frame(prodPts(object))
        pf=model.frame(prodFn( object))
            
        msy=pf[which.min(abs(pf$yield-max(pf$yield))), ]
        pf2=subset(pf,ssb<=msy$ssb)
            
        if (rm.neg) pt=subset(pt,yield>0)
        yield=median(pt$yield,na.rm=TRUE)
        
        # Use approx for linear interpolation
        hat=approx(x=pf2$yield,y=pf2$ssb,xout=yield)$y
            
        rtn=data.frame(ssb=hat,yield=yield)
        return(rtn)})

