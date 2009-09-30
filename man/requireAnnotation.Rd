\name{requireAnnotation}
\Rdversion{1.1}
\alias{requireAnnotation}
\title{
  Helper function to load packages.
}
\description{
  This function checkes the existence of a given package and loads it if
  available. If the package is not available, the function checks its
  availability on BioConductor, downloads it and installs it.
}
\usage{
requireAnnotation(pkgname, lib=.libPaths()[1], verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pkgname}{
    character. Package name (usually an annotation package).
  }
  \item{lib}{
    character. Path where to install packages at.
    }
  \item{verbose}{
    logical. Verbosity flag.
}
}
\value{
  Logical: TRUE if package is available or FALSE if package unavailable
  for download.
}
\author{
  Benilton Carvalho
}

\seealso{
  install.packages
}
\examples{
\dontrun{
requirePackage("pd.mapping50k.xba240")
}
}

\keyword{data}

