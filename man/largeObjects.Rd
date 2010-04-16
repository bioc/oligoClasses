\name{initializeBigMatrix}
\alias{initializeBigVector}
\alias{initializeBigMatrix}

\title{
  Initialize big matrices/vectors.
}
\description{
  Initialize big matrices or vectors appropriately (conditioned on the
  status of support for large datasets - see Details).
}
\usage{
initializeBigMatrix(name, nr, nc, vmode = "integer", initdata = NA)
initializeBigVector(name, n, vmode = "integer", initdata = NA)
}

\arguments{
  \item{name}{prefix to be used for file stored on disk}
  \item{nr}{number of rows}
  \item{nc}{number of columns}
  \item{n}{length of the vector}
  \item{vmode}{mode - "integer", "double"}
  \item{initdata}{Default is NA} %%FIXME
}
\details{

  These functions are meant to be used by developers. They provide means
  to appropriately create big vectors or matrices for packages like
  oligo and crlmm (and friends). These objects are created conditioned
  on the status of support for large datasets.
  
}
\value{
  If the 'ff' package is loaded (in the search path), then an 'ff'
  object is returned. A regular R vector/matrix is returned otherwise.
}

\examples{
x <- initializeBigVector("test", 10)
class(x)
x
if (isPackageLoaded("ff"))
  finalizer(x) <- "delete"
rm(x)
}

\keyword{ manip }

