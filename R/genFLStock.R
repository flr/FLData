###############################################################################
# EJ&CM(2012821)
# Methods to add uncertainty to FLQuant objects
# NOTE #1:
# NOTE #2:
###############################################################################

#==============================================================================
# 
#==============================================================================

#' Methods to genetate FLStock objects
#'
#' Some additional details about this S4 generic and its methods.
#' The extra blank line between this section and the title is
#' critical for roxygen2 to differentiate the title from the
#' description section.
#'
#' @param object an FLQuant
#' @param rec an FLQuant
#' @param catch.n an FLQuant
#' @param harvest an FLQuant
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
#' @rdname au-methods
#'
#' @examples
#' data(ple4)
#' sim.F <- genFLQuant(harvest(ple4))
setGeneric("au", function(object, R, C, F, ...) standardGeneric("au"))

#' @rdname au-methods
#' @aliases au,FLStock,FLQuant,FLQuant,missing-method
setMethod("au", c("FLStock", "FLQuant", "FLQuant", "missing"), function(object, R, C, F, ...){
	cat("Not implemented yet\n")
})

#' @rdname au-methods
#' @aliases au,FLStock,missing,FLQuant,FLQuant-method
setMethod("au", c("FLStock", "missing", "FLQuant", "FLQuant"), function(object, R, C, F, ...){
	cat("Not implemented yet\n")
})

#' @rdname au-methods
#' @aliases au,FLStock,FLQuant,missing,FLQuant-method
setMethod("au", c("FLStock", "FLQuant", "missing", "FLQuant"), function(object, R, C, F, ...){
	# requires checking dimensions
	if(!identical(dim(catch.n(object))[-c(1,6)], dim(R)[-c(1,6)])) stop("Recruitment vector must have consistent dimensions with the stock object")
	if(!identical(dim(catch.n(object))[-6]    , dim(F)[-6])) stop("Harvest matrix must have consistent dimensions with the stock object")
	if(dim(R)[6]!=dim(R)[6]) stop("R and F must have the same number of iterations")

	# get dims and set flq
	dms <- dims(object)
	nages <- dms$age
	nyrs <- dms$year
	niters <- dim(R)[6]
	flq <- FLQuant(dimnames=dimnames(F))
	
	# compute cumulative Z
	Z <- FLCohort(F + m(object))
	Z[] <- apply(Z, c(2:6), cumsum)

	# expand variability into [N] by R*[F] 
	Ns <- FLCohort(R[rep(1,nages)])
	Ns <- Ns*exp(-Z)
	Ns <- as(Ns, "FLQuant")

	# Update object
	stock.n(object) <- flq
	# [R]
	stock.n(object)[1] <- R
	# [N]
	stock.n(object)[-1,-1] <- Ns[-nages,-nyrs] 
	# plus group
	stock.n(object)[nages,-1] <- Ns[nages,-1] + stock.n(object)[nages,-1]
	stock(object) <- computeStock(object)
	# [F]
	harvest(object) <- F
	# [C]
	Z <- harvest(object) + m(object)
	Cs <- harvest(object)/Z*(1-exp(-Z))*stock.n(object) 
	catch.n(object) <- Cs
	catch(object) <- computeCatch(object)
	# [L] & [D] rebuilt from C
	# Ds=D/(D+L)*Cs where Cs is the simulated catch
	D <- discards.n(object)
	L <- landings.n(object)
	discards.n(object) <- D/(D+L)*Cs
	discards(object) <- computeDiscards(object)
	landings.n(object) <- Cs - discards.n(object)
	landings(object) <- computeLandings(object)
	# out
	object
})


