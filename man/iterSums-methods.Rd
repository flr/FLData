\docType{methods}
\name{iterSums}
\alias{iterSums}
\alias{iterSums,FLQuant-method}
\title{Calculate the sums accross iterations}
\arguments{
  \item{object}{an FLQuant with iters}

  \item{...}{Additional argument list that might not ever
  be used.}
}
\value{
  an FLQuant
}
\description{
  Calculate the sums accross iterations
}
\examples{
data(ple4)
genFLQuant(harvest(ple4), method = "ac")
}
\seealso{
  \code{\link{print}} and \code{\link{cat}}
}

