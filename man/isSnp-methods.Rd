\name{isSnp-methods}
\docType{methods}
\alias{isSnp-methods}
\alias{isSnp}
\alias{isSnp,character-method}
\alias{isSnp,SnpSet-method}
\alias{isSnp,gSet-method}
\alias{isSnp,GenomeAnnotatedDataFrame-method}
\title{ Methods for Function isSnp in package oligoClasses~~}
\description{
 ~~ Methods for function \code{isSnp} in package \pkg{oligoClasses} ~~
}
\section{Methods}{
\describe{

  Return an indicator for whether the marker is polymorphic (value 1) or
  nonpolymorphic (value 0).

\item{\code{signature(object = "character", pkgname = "character")}}{
  Return an indicator for whether the vector of marker identifiers in
  \code{object} is polymorphic.  \code{pkgname} must be one of the
  supported annotation packages specific to the platform.
}

\item{\code{signature(object = "eSet", pkgname = "ANY")}}{
  If 'isSnp' is included in \code{fvarLabels(object)}, an indicator for
  polymorphic markers is returned. Otherwise, an error is thrown.
}

\item{\code{signature(object = "GenomeAnnotatedDataFrame", pkgname =
  "ANY")}}{ Accessor for indicator of whether the marker is polymorphic.
  If annotation was not available due to a missing or non-existent
  annotation package, the value returned by the accessor will be a
  vector of zero's.

}
}
}

\keyword{methods}

