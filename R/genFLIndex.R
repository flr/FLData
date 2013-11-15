###############################################################################
# EJ&CM(2012821)
# Methods to create index FLQuants from stock.n
# NOTE #1:
# NOTE #2:
###############################################################################

#==============================================================================
# 
#==============================================================================

#' Methods to create index FLQuants from stock.n
#'
#' Some additional details about this S4 generic and its methods.
#' The extra blank line between this section and the title is
#' critical for roxygen2 to differentiate the title from the
#' description section.
#'
#' @param object an FLQuant containing stock numbers at age
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant?
#' 
#' @seealso \code{\link{stock.n}} and \code{\link{FLQuant}}
#' 
#' @export
#' @docType methods
#' @rdname genFLIndex-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("genFLIndex", function(object, ...) standardGeneric("genFLIndex"))


#' @rdname genFLIndex-methods
#' @aliases genFLIndex,FLQuant-method
setMethod("genFLIndex", c("FLQuant"), 
    function(object, cv = 0.2, 
                     catchability = "flat", 
                     niter = 250) {
      # use log transform, to be expanded on later versions
      mu <- log(object)
      
      if(method == "ac") {
        Rho <- cor(t(mu[drop = TRUE]))
        flq <- mvrnorm(niter * dim(mu)[2], rep(0, nrow(Rho)), cv^2 * Rho)
        mu <- propagate(mu, niter)
        flq <- FLQuant(c(t(flq)), dimnames = dimnames(mu))
        flq <- exp(mu + flq)
      }

      
      return(flq)
})

