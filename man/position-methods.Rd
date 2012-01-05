\name{position-methods}
\docType{methods}
\alias{position-methods}
\alias{position}
\alias{position,AnnotatedDataFrame-method}
\alias{position,SnpSet-method}
\alias{position,gSet-method}
\alias{position,GenomeAnnotatedDataFrame-method}
\title{Methods for function position in Package oligoClasses}
\description{
   Methods for function \code{position} in package \pkg{oligoClasses}
}
\section{Methods}{

  The methods for \code{position} extracts the physical position stored
  as an integer for each marker in a \code{eSet}-derived class or a
  \code{AnnotatedDataFrame}-derived class.

\describe{

\item{\code{signature(object = "AnnotatedDataFrame")}}{
  Accessor for physical position.
}

\item{\code{signature(object = "eSet")}}{
  If 'position' is included in \code{fvarLabels(object)}, the physical
  position will be returned. Otherwise, an error is thrown.
}

\item{\code{signature(object = "GenomeAnnotatedDataFrame")}}{
  Accessor for physical position.  If annotation was not available  due
  to a missing or non-existent annotation package, the value  returned
  by the accessor will be a vector of zero's.
}
}
}
\keyword{methods}

