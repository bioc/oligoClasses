\name{eSet-methods}
%\docType{methods}
\alias{chromosome}
\alias{chromosome,eSet-method}
\alias{chromosome}
\alias{chromosome<-}
\alias{chromosome<-,eSet-method}
\alias{chromosome,AnnotatedDataFrame-method}
\alias{position}
\alias{position,eSet-method}
\alias{position,AnnotatedDataFrame-method}
\alias{isSnp}
\alias{isSnp,eSet-method}
\title{Accessors for eSet extensions}
\description{

	Accessors for variables stored in the featureData slot of a
	class inheriting from eSet.

      }

\usage{
chromosome(object)
chromosome(object) <- value
position(object)
isSnp(object)
}

\arguments{
  \item{object}{A \code{eSet} object or \code{AnnotatedDataFrame}.}
  \item{value}{A \code{integer} vector.}
}
\keyword{methods}
\keyword{manip}
