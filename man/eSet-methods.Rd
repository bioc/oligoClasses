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
\alias{isSnp,eSet,ANY-method}
\alias{isSnp,character,character-method}
\title{Accessors for eSet extensions}
\description{

	Accessors for variables stored in the featureData slot of a
	class inheriting from eSet.

      }

\usage{
chromosome(object,na.rm=FALSE)
chromosome(object) <- value
position(object,na.rm=FALSE)
isSnp(object,pkgname)
##snpNames(object)
}

\arguments{
  \item{object}{A \code{eSet} object or \code{AnnotatedDataFrame}.}
  \item{na.rm}{Logical.  Whether to exlude NA's.}
  \item{value}{A \code{integer} vector.}
  \item{pkgname}{A \code{character} string indicating the annotation package.}
}
\value{
  \code{position} returns an integer vector of the genomic position
  position (units are base pairs).

  \code{chromosome} returns an integer vector of the chromosome.  See
  \code{chromosome2integer} for the integer coding of non-autosomal
  chromosomes.

  \code{isSnp} returns a logical vector that indicates which markers are
  polymorphic. If \code{object} is character vector of feature names,
  the \code{pkgname} must be supplied. For nonpolymorphic markers,
  \code{isSnp} evaluates to FALSE.

}
\seealso{\code{\link{chromosome2integer}}, \code{\link{annotationPackages}}}

\keyword{methods}
\keyword{manip}
