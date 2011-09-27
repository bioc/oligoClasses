\name{CopyNumberSet-methods}
\docType{methods}
\alias{CopyNumberSet-methods}
\alias{copyNumber}
\alias{copyNumber<-}
\alias{cnConfidence}
\alias{cnConfidence<-}
\title{Methods for class CopyNumberSet.}
\description{
	Accessors and CopyNumberSet
}

\usage{

copyNumber(object, ...)
cnConfidence(object)
copyNumber(object) <- value
cnConfidence(object) <- value

}

\arguments{

\item{object}{\code{CopyNumberSet} object or derived class}
\item{\dots}{Ignored for \code{CopyNumberSet} and \code{oligoSnpSet}. }
\item{value}{matrix}
}

\value{

	\code{copyNumber} returns a matrix of copy number estimates.
	\code{cnConfidence} returns a matrix of confidence scores for
	the copy number estimates.

}
\keyword{methods}
\keyword{manip}

