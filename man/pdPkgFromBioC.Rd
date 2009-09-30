\name{pdPkgFromBioC}
\Rdversion{1.1}
\alias{pdPkgFromBioC}

\title{
  Get packages from BioConductor.
}
\description{
  This function checks if a given package is available on BioConductor
  and installs it, in case it is.
}
\usage{
pdPkgFromBioC(pkgname, lib = .libPaths()[1], verbose = TRUE)
}

\arguments{
  \item{pkgname}{
    character. Name of the package to be installed.
}
  \item{lib}{
    character. Path where to install the package at.
}
  \item{verbose}{
    logical. Verbosity flag.
}
}
\details{
  Internet connection required.
}
\value{
  Logical: TRUE if package was found, downloaded and installed; FALSE otherwise.
}
\author{
  Benilton Carvalho
}

\seealso{
  download.packages
}
\examples{
\dontrun{
pdPkgFromBioC("pd.mapping50k.xba240")
}
}
\keyword{data}

