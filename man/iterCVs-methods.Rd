\docType{methods}
\name{iterCVs}
\alias{iterCVs}
\alias{iterCVs,FLQuant-method}
\title{Calculate the CVs accross iterations}
\arguments{
  \item{object}{an FLQuant with iters}

  \item{...}{Additional argument list that might not ever
  be used.}
}
\value{
  an FLQuant
}
\description{
  Calculate the CVs accross iterations
}
\examples{
data(ple4)
genFLQuant(harvest(ple4), method = "ac")
}
\seealso{
  \code{\link{print}} and \code{\link{cat}}
}

