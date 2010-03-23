\name{isPackageLoaded}
\alias{isPackageLoaded}

\title{
  Check if package is loaded.
}
\description{
  Checks if package is loaded.
}
\usage{
isPackageLoaded(pkg)
}

\arguments{
  \item{pkg}{
    Package to be checked.
}
}
\details{
  Checks if package name is in the search path.
}
\value{
  Logical.
}

\seealso{
  search
}
\examples{
isPackageLoaded("oligoClasses")
isPackageLoaded("ff")
isPackageLoaded("snow")
}

\keyword{manip}

