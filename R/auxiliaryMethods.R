#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname iterMedians-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("iterMedians", function(object, ...) standardGeneric("iterMedians"))

#' @rdname iterMedians-methods
#' @aliases iterMedians,FLQuant-method
setMethod("iterMedians", "FLQuant", function(object, ...){
	return(apply(object, c(1:5), median, na.rm = FALSE))
})

#' Calculate the sums accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname iterSums-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("iterSums", function(object, ...) standardGeneric("iterSums"))

#' @rdname iterSums-methods
#' @aliases iterSums,FLQuant-method
setMethod("iterSums", "FLQuant", function(object, ...){
	return(apply(object, c(1:5), sum, na.rm = FALSE))
})

#' Calculate the CVs accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname iterCVs-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("iterCVs", function(object, ...) standardGeneric("iterCVs"))

#' @rdname iterCVs-methods
#' @aliases iterCVs,FLQuant-method
setMethod("iterCVs", "FLQuant", function(object, ...){
	return(sqrt(iterVars(object))/iterMeans(object))
})


