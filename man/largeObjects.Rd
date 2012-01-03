\name{initializeBigMatrix}
\alias{initializeBigVector}
\alias{initializeBigMatrix}
\alias{initializeBigArray}

\title{
  Initialize big matrices/vectors.
}
\description{
  Initialize big matrices or vectors appropriately (conditioned on the
  status of support for large datasets - see Details).
}
\usage{
initializeBigMatrix(name=basename(tempfile()), nr=0L, nc=0L, vmode = "integer", initdata = NA)
initializeBigVector(name=basename(tempfile()), n=0L, vmode = "integer",
  initdata = NA)
initializeBigArray(name=basename(tempfile()), dim=c(0L,0L,0L),
  vmode="integer", initdata=NA)
}

\arguments{
  \item{name}{prefix to be used for file stored on disk}
  \item{nr}{number of rows}
  \item{nc}{number of columns}
  \item{n}{length of the vector}
  \item{vmode}{mode - "integer", "double"}
  \item{initdata}{Default is NA} %%FIXME
  \item{dim}{Integer vector indicating the dimensions of the array to initialize}
}
\details{

  These functions are meant to be used by developers. They provide means
  to appropriately create big vectors or matrices for packages like
  oligo and crlmm (and friends). These objects are created conditioned
  on the status of support for large datasets.

}
\value{
  If the 'ff' package is loaded (in the search path), then an 'ff'
  object is returned. A regular R vector or array is returned otherwise.
}

\examples{
x <- initializeBigVector("test", 10)
class(x)
x
if (isPackageLoaded("ff"))
  finalizer(x) <- "delete"
rm(x)
initializeBigMatrix(nr=5L, nc=5L)
initializeBigArray(dim=c(10, 5, 3))
}

\keyword{ manip }

