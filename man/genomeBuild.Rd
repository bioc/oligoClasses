\name{genomeBuild}
\alias{genomeBuild}
\alias{genomeBuild<-}
\alias{genomeBuild,DBPDInfo-method}
\alias{genomeBuild,FeatureSet-method}

\title{Genome Build Information}
\description{
  Returns the genome build. This information comes from the
  annotation package and is given as an argument during the package
  creation process.
}
\usage{
genomeBuild(object)
}

\value{ character string  }

\note{Supported builds are UCSC genome builds are 'hg18' and 'hg19'.}

\arguments{
  \item{object}{Supported objects include \code{PDInfo},
  \code{FeatureSet}, and any \code{gSet}-derived or
  \code{eSetList}-derived object.  }
}

\examples{
  showMethods("genomeBuild", where="package:oligoClasses")
}

\keyword{manip}

