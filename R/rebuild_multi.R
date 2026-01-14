# =============================================================================
# Rebuild Methods for Multiple Scenarios
# =============================================================================
# Note: Generic definitions are in generic.R
# Note: Main rebuild methods are in rebuild_methods.R

# =============================================================================
# rebuild2 Function
# =============================================================================

#' Rebuild trajectories for multiple scenarios (version 2)
#' 
#' @description Projects rebuilding trajectories from different initial SSB levels using version 2 algorithm
#' 
#' @param object An FLBRP object
#' @param targetF Target fishing mortality (default is 0)
#' @param targetSSB Target SSB (default is MSY SSB)
#' @param initial Number of initial depletion levels to evaluate (default 100)
#' @param minVal Minimum depletion level (default 1e-6)
#' @param maxVal Maximum depletion level (default 1)
#' @param growthRate Growth rate parameter for SSB sequence (default 0.25)
#' @return FLStock object with rebuilding trajectories
#' @export
rebuild2<-function(object, 
                   targetF  =computeRefpts(object)["msy","harvest"]*0,
                   targetSSB=computeRefpts(object)["msy","ssb"],
                   initial=100, minVal=1e-6, maxVal=1,growthRate = 0.25) {
            
            # Input validation
            if (!is(object, "FLBRP")) 
              stop("object must be an FLBRP object")
            
            if (!(is.numeric(initial) && initial > 0))
              stop("initial must be positive integer")
            
            nits=max(as.integer(initial),1)
            
            if (!all(sapply(list(growthRate, minVal, maxVal), is.numeric))) 
              stop("growthRate, minVal, and maxVal must be numeric")
            
            if (minVal > maxVal)
              stop("minVal must be less than maxVal")
            
            if (nits==1) minVal=maxVal=initial
            
            # Generate F by SSB
            rfs=propagate(refpts(object)["msy",],nits)
            dimnames(rfs)$refpt[1]="ssb"
            rfs["ssb","ssb"]=rfs["ssb","ssb"]*seq(minVal^growthRate, maxVal^growthRate, 
                                                         length.out=nits)^(1/growthRate)
            rfs["ssb",-4]=NA

            eql        =object
            refpts(eql)=rfs
            refpts(eql)=computeRefpts(eql)
            
            fbar(eql)=FLQuant(1,dimnames=list(year=seq(101+dim(m(eql))[1]),iter=nits))
            fbar(eql)=fbar(eql)%*%refpts(eql)["ssb","harvest"]
            
            eql=brp(eql)
                    
            # Project stock
            stk =as(eql, "FLStock")
            fbar(eql)[,-seq(dim(m(eql))[1])] = targetF
            
            stk =ffwd2(stk, f=fbar(eql)[,-1], sr=eql)
            stk =stk[,-seq(dim(m(eql))[1])]
            
            # Rename years
            stk = qapply(stk, function(x) {
              dimnames(x)$year = seq(length(dimnames(x)$year))
              x})
            
            return(stk)}
