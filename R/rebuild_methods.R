# =============================================================================
# Rebuild Methods for FLBRP and FLStock
# =============================================================================
# Note: Generic definitions are in generic.R

# =============================================================================
# Rebuild Methods for FLBRP and FLStock
# =============================================================================
# Note: Helper functions (interp, tryIt) are in helpers.R
# Generic definitions are in generic.R

# =============================================================================
# rebuild Methods
# =============================================================================

#' @importFrom FLasher fwd
#' @rdname rebuild
#' @export
setMethod("rebuild", signature(object = "FLBRP"),
  function(object, 
           targetF   = refpts(object)["msy","harvest"] * 0,
           targetSSB = refpts(object)["msy","ssb"],
           initial = 100, growthRate = 0.25, minVal = 1e-6, maxVal = 1,
           burnin = 20, truncate = TRUE) {
    
    # Input validation
    if (!is(object, "FLBRP")) {
      stop("object must be an FLBRP object")
    }
    
    if (!all(sapply(list(initial, burnin), function(x) is.numeric(x) && x > 0))) {
      stop("initial and burnin must be positive integers")
    }
    nits=max(as.integer(initial),1)
    
    if (!all(sapply(list(growthRate, minVal, maxVal), is.numeric))) {
      stop("growthRate, minVal, and maxVal must be numeric")
    }
    
    if (minVal > maxVal) {
      stop("minVal must be less than maxVal")
    }
    if (nits==1) minVal=maxVal=initial
    
    # Generate SSB sequence
    targetSSB = c(targetSSB) * seq(minVal^growthRate, maxVal^growthRate, 
                                  length.out = nits)^(1/growthRate)
    targetF   = c(targetF)
    
    # Setup equilibrium
    eql = object
    fbar(eql)[] = 0.2
    
    # Create target biomass array
    btar = FLQuant(rep(targetSSB, each = dim(fbar(eql))[2]),
                   dimnames = dimnames(propagate(ssb(eql), nits)))
  
    # Project stock
    stk = propagate(as(eql, "FLStock"), nits)
    stk = FLasher::fwd(stk, ssb_end = btar[,-seq(dims(stk)[["min"]]+2)], sr = eql)
   
    # Apply target F
    ftar = fbar(stk) %=% targetF

    #stk1 = window(fwd(iter(stk,2), f=ftar[,-seq(burnin),,,,1], sr=eql),start=burnin+1)
    #rcvYr=min(as.numeric(dimnames(ssb(stk1))$year[c(ssb(stk1)>max(targetSSB))]))
    #stk  =window(stk,end=rcvYr)

    stk = ffwd2(stk, f = ftar[,-seq(burnin)], sr = eql)
  
    # Post-process
    if (truncate) {
      stk = stk[,-seq(burnin)]
    }
    
    # Rename years
    stk = qapply(stk, function(x) {
      dimnames(x)$year = seq(length(dimnames(x)$year))
      x
    })
    
    return(stk)
  })

#' @rdname rebuild
#' @export
setMethod("rebuild", signature(object = "biodyn"),
  function(object, target = refpts(object)["bmsy"], initial = 100, 
           growthRate = 0.3, minVal = 1e-6, maxVal = 1, burnin = 20, truncate = TRUE) {
    
    # Input validation
    if (!is(object, "biodyn")) {
      stop("object must be a biodyn object")
    }
    
    if (!is.numeric(target) || length(target) != 1) {
      stop("target must be a single numeric value")
    }
    
    if (!all(sapply(list(initial, burnin), function(x) is.numeric(x) && x > 0))) {
      stop("initial and burnin must be positive integers")
    }
    
    if (!all(sapply(list(growthRate, minVal, maxVal), is.numeric))) {
      stop("growthRate, minVal, and maxVal must be numeric")
    }
    
    if (minVal >= maxVal) {
      stop("minVal must be less than maxVal")
    }
    
    # Get BMSY
    bmsy <- c(refpts(object)["bmsy"])
    
    # Create propagated object
    rtn <- propagate(object, initial)
    
    # Create stock projection
    target_seq <- c(target) * seq(minVal^growthRate, maxVal^growthRate, 
                                length.out = initial)^(1/growthRate)
    rtn@stock <- FLQuant(rep(target_seq, each = dim(rtn)[2]), 
                       dimnames = dimnames(stock(rtn)))
    rtn <- fwd(rtn, catch = catch(rtn)[,-1] %=% 0.0)
    
    # Post-process
    if (truncate) {
      rtn <- rtn[,-seq(burnin)]
    }
    
    # Rename years
    rtn <- qapply(rtn, function(x) {
      dimnames(x)$year <- seq(length(dimnames(x)$year))
      x
    })
    
    return(rtn)
  })

# =============================================================================
# rebuildTime Methods
# =============================================================================

#' @rdname rebuildTime
#' @export
setMethod("rebuildTime", signature(object = "FLStock"),
  function(object) {

    # Extract SSB data
    df <- as.data.frame(ssb(object), drop = TRUE)
    iters <- sort(an(unique(df$iter)))

    # Calculate BMSY and scale SSB
    bmsy <- c(ssb(object)[,1,,,,dim(object)[6]])
    df$biomass <- df$data/bmsy
    df$initial <- rep(c(ssb(object[,1]))[iters]/bmsy, each=dim(object)[2])
    df <- na.omit(df)

    # Interpolate results
    return(interp(df))})

#' @rdname rebuildTime
#' @export
setMethod("rebuildTime", signature(object = "biodyn"), 
  function(object, target = refpts(object)["bmsy"], initial = 100, 
           growthRate = 0.3, minVal = 1e-6, maxVal = 1, nx = 101) {

    # Input validation
    if (!is(object, "biodyn")) {
      stop("object must be a biodyn object")
    }
    
    if (!is.numeric(target) || length(target) != 1) {
      stop("target must be a single numeric value")
    }
    
    if (!all(sapply(list(initial, nx), function(x) is.numeric(x) && x > 0))) {
      stop("initial and nx must be positive integers")
    }
    
    if (!all(sapply(list(growthRate, minVal, maxVal), is.numeric))) {
      stop("growthRate, minVal, and maxVal must be numeric")
    }
    
    if (minVal >= maxVal) {
      stop("minVal must be less than maxVal")
    }
    
    # Get BMSY
    bmsy = c(refpts(object)["bmsy"])
    
    # Create propagated object
    rtn = propagate(object, initial)
    
    # Create stock projection
    target_seq = c(target) * seq(minVal^growthRate, maxVal^growthRate, 
                                length.out = initial)^(1/growthRate)
    rtn@stock = FLQuant(rep(target_seq, each = dim(rtn)[2]), 
                       dimnames = dimnames(stock(rtn)))
    rtn = fwd(rtn, catch=catch(rtn)[,-1] %=% 0.0)
    
    # Transform data
    df = as.data.frame(stock(rtn), drop = TRUE)
    df$initial = c(stock(rtn)[,1])[an(df$iter)]
    df = df[,-2]
    names(df) = c("year","biomass","initial")
    df[,-1]=df[,-1]/bmsy

    # Interpolate results
    return(interp(df))
  })

rebuildTimeFLBRP<-function(object, 
                   targetF  =computeRefpts(object)["msy","harvest"]*0,
                   targetSSB=computeRefpts(object)["msy","ssb"],
                   minVal   = 0.05) {
  
  # Input validation
  if (!is(object, "FLBRP")) 
    stop("object must be an FLBRP object")
  
  if (minVal < 0) 
    stop("minVal must be greater than 0")
  if (minVal > 1) 
    stop("maxVal must be less than 1")
  
  finitial=computeRefpts(object)["crash","harvest"]*0.99
  if (!is.na(c(finitial))){
        refpts(object)["crash","harvest"]*0.99}else{
          refpts(object)=refpts(object)[1,]
          dimnames(refpts(object))[[1]]="ssb"
          refpts(object)[]=NA
          refpts(object)[1,1]=targetSSB*1e-4
          finitial=computeRefpts(object)[1,"harvest"]
        }
  
  fbar(object)[]=finitial
  object        =brp(object)
  stk           =as(object,"FLStock")
  
  if (dims(stk)["min"]==0)
     rtn=tryIt(FLCore::ffwd(stk, fbar=fbar(object)[,-1]%=%0, sr=object))
  else
    rtn=tryIt(ffwd(stk, f=fbar(object)[,-seq(dims(stk)[["min"]]+1)]%=%0, sr=object))
    #rtn=tryIt(fwd(stk, f=fbar(object)[,-seq(dims(stk)[["min"]]+1)]%=%0, sr=object))
  
  df=as.data.frame(ssb(rtn[,-seq(dims(stk)[["min"]]+1)]) %/% refpts(object)["msy","ssb"], drop = TRUE)[-1,]
  df=subset(df, data <= 1)
  df$year = max(df$year)-df$year+1
  names(df)=c("T_recover","initial")
  df=rbind(df[order(df$initial),],data.frame(initial=1,T_recover=0))
  return(df)}


