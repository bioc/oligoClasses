\name{CopyNumberSet-methods}
\docType{methods}
\alias{CopyNumberSet-methods}
\alias{copyNumber}
<<<<<<< HEAD
\alias{copyNumber<-}
\alias{cnConfidence}
\alias{cnConfidence<-}
=======
%\alias{copyNumber,CopyNumberSet-method}
\alias{copyNumber<-}
%\alias{copyNumber<-,CopyNumberSet,matrix-method}
\alias{cnConfidence}
\alias{cnConfidence<-}
%\alias{cnConfidence,CopyNumberSet-method}
%\alias{cnConfidence<-,CopyNumberSet,matrix-method}
\alias{initialize,CopyNumberSet-method}
>>>>>>> Export CopyNumberSet class
\title{  Methods for class CopyNumberSet.}
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

