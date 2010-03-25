\name{is.ffmatrix}
\alias{is.ffmatrix}

\title{
Check if object is an ff-matrix object.
}
\description{
Check if object is an ff-matrix object.
}
\usage{
is.ffmatrix(object)
}

\arguments{
  \item{object}{object to be checked}
}

\value{
  Logical.
}

\note{
  This function is meant to be used by developers.
}

\examples{
if (isPackageLoaded("ff")){
  x1 <- ff(vmode="double", dim=c(10, 2))
  is.ffmatrix(x1)
}
x1 <- matrix(0, nr=10, nc=2)
is.ffmatrix(x1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
