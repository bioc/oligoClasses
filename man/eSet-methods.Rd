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
\alias{snpNames}
\alias{snpNames,eSet-method}
\title{Accessors for eSet extensions}
\description{

	Accessors for variables stored in the featureData slot of a
	class inheriting from eSet.

      }

\usage{
chromosome(object,na.rm=FALSE)
chromosome(object) <- value
position(object,na.rm=FALSE)
isSnp(object)
snpNames(object)
}

\arguments{
  \item{object}{A \code{eSet} object or \code{AnnotatedDataFrame}.}
  \item{na.rm}{Logical.  Whether to exlude NA's.}
  \item{value}{A \code{integer} vector.}
}
\keyword{methods}
\keyword{manip}
