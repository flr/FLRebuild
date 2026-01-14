#' @importFrom FLCore units
setGeneric("expand", useAsDefault=Matrix::expand)

setMethod('expand', signature(x='FLArray'),
        function(x, ..., fill=TRUE) {

            args <- list(...)
            dnx <- dimnames(x)

            # dimension names
            nargs <- names(args)
            qnames <- names(dnx)

            # check input names match dimnames
            if(!all(nargs %in% qnames))
              stop(paste("Wrong dimension name provided: ", nargs[!nargs%in%qnames]))

            # turn into characters
            select <- lapply(args, as.character)

            # match specified dimensions and dimnames
            dimnames <- dnx

            # new dimnames
            dimnames[names(select)] <- select

            # output object
            res <- new(class(x), array(as.numeric(NA), dimnames=dimnames,
                                       dim=unlist(lapply(dimnames, length))), units=units(x))

            # IF !fill, return
            if(!fill)
              return(res)

            # ANY dim with new names?
            new <- unlist(Map(function(x,y) !all(x %in% y), x=dnx, y=dimnames))

            # ASSIGN new dimnames in 'new' dims
            dnx[new] <- dimnames(res)[new]

            # list names to match '[<-' signature
            names(dnx) <- c('i', 'j', 'k', 'l', 'm', 'n')

            return(do.call('[<-', c(list(x=res, value=x), dnx)))
          })


#' Forward Project Stock with Age Checking
#'
#' Wrapper that ensures recruitment age compatibility, runs ffwd, then trims ages if needed.
#' Automatically expands to age 0 and resets result to original age range.
#'
#' @param stk FLStock object
#' @param control FLPar or control object for ffwd
#' @param ... Additional arguments passed to ffwd
#'
#' @importFrom FLCore ffwd
#' @export
#'
#' @examples
#' \dontrun{
#' library(FLCore)
#' # Create an FLStock with first age > 0
#' data(ple4)
#' stk <- ple4
#' # Create control FLPar
#' control <- FLPar(fbar=0.2)
#' # Forward project
#' out <- ffwd2(stk, control)
#' }

ffwd2 <- function(object, control, ...) {
  orig = dimnames(object)$age
  minAge = dims(object)[["min"]]

  if (minAge > 0) {
    object = suppressWarnings(FLCore::expand(object, age = seq(0, dims(object)[["max"]])))

    object = suppressWarnings(qapply(object, function(x) {x[is.na(x)] = 0; x}))
    
    for (iAge in minAge:1)
      stock.n(object)[iAge, -dim(object)[2]] = stock.n(object)[iAge+1, -1]
  }
 
  res = suppressWarnings(FLCore::ffwd(object, control, ...))
  
  if (minAge > 0)
    res = trim(res, age = orig)
  
  return(res)}
