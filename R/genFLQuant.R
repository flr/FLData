###############################################################################
# EJ&CM(2012821)
# Methods to add uncertainty to FLQuant objects
# NOTE #1:
# NOTE #2:
###############################################################################

#==============================================================================
# 
#==============================================================================

#' Methods to add uncertainty to FLQuant objects
#'
#' Some additional details about this S4 generic and its methods.
#' The extra blank line between this section and the title is
#' critical for roxygen2 to differentiate the title from the
#' description section.
#'
#' @param object an FLQuant
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
#' @rdname getAcor-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("getAcor", function(object, ...) standardGeneric("getAcor"))

#' @rdname getAcor-methods
#' @aliases getAcor,FLQuant-method
setMethod("getAcor", c("FLQuant"), function(object, tf=log, ...) {
		mu <- log(object)
		Rho <- cor(t(mu[drop = TRUE]))
		return(Rho)
})

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
#' @rdname genFLQuant-methods
#'
#' @examples
#' data(ple4)
#' sim.F <- genFLQuant(harvest(ple4), method = "ac")
setGeneric("genFLQuant", function(object, ...) standardGeneric("genFLQuant"))


#' @rdname genFLQuant-methods
#' @aliases genFLQuant,FLQuant-method
setMethod("genFLQuant", c("FLQuant"), function(object, cv = 0.2, method = "rw", niter = 250) {
  # use log transform, to be expanded on later versions
	mu <- log(object)
	if(method == "ac") {
		Rho <- cor(t(mu[drop = TRUE]))
		flq <- mvrnorm(niter * dim(mu)[2], rep(0, nrow(Rho)), cv^2 * Rho)
		mu <- propagate(mu, niter)
		flq <- FLQuant(c(t(flq)), dimnames = dimnames(mu))
		flq <- exp(mu + flq)
	}
	if(method == "rw") {
		n.ages  <- dim(mu)[1]
	  	n.years <- dim(mu)[2]
		n <- n.ages * n.years
		# set up lcs to extract posterior means
		B = diag(n)
		B[B==0] <- NA
		lcs = inla.make.lincombs(Predictor = B)
		# treat mu as a GMRF model - 
		# an independant random walk for each age (same variance accross ages)
		form <- x ~ f(years, model = 'rw1', replicate = ages)
		data <- list(x = c(mu), years = rep(1:n.years, each = n.ages), ages  = rep(1:n.ages, n.years))
		result <- inla(form, data = data, control.predictor = list(compute = TRUE), 
                       lincomb = lcs, control.inla = list(lincomb.derived.correlation.matrix = TRUE))
		# the covariance of the fitted RW
		RW.cov <- result $ misc $ lincomb.derived.correlation.matrix
		# two options for the mean:
		#  1) use the mean estimate of RW process from INLA
		#     - this is potentially very smooth and lacking in strucure
		#	mu.hat <- result $ summary.linear.predictor $ mean
		#	flq <- mvrnorm(niter, mu.hat, cv^2 * RW.cov)
		#  2) use the original data and add the noise to that
		#  2 is more consistent with ac method and always maintains data structure
		flq <- exp(mvrnorm(niter, c(mu), cv^2 * RW.cov))
		flq <- FLQuant(c(t(flq)), dimnames = dimnames(propagate(mu, niter)))
	}
	units(flq) <- units(object)
	return(flq)
})

