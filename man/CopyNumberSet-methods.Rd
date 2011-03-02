\name{CopyNumberSet-methods}
\docType{methods}
\alias{CopyNumberSet-methods}
\alias{copyNumber, ...}
\alias{copyNumber,CopyNumberSet-method}
\alias{copyNumber<-}
\alias{copyNumber<-,CopyNumberSet,matrix-method}
\alias{cnConfidence}
\alias{cnConfidence<-}
\alias{cnConfidence,CopyNumberSet-method}
\alias{cnConfidence<-,CopyNumberSet,matrix-method}
\alias{initialize,CopyNumberSet-method}
\title{  Methods for class CopyNumberSet.}
\description{
	Accessors and CopyNumberSet
}

\usage{

copyNumber(object)
cnConfidence(object)

}

\arguments{

\item{object}{\code{CopyNumberSet} object}
\item{\dots}{Ignored for CopyNumberSet and oligoSnpSet. For CNSet,
specification of the marker indices and/or column indices is required
through arguments \code{i} and \code{j}, respectively.
}
}

\value{

	\code{copyNumber} returns a matrix of copy number estimates.
	
	\code{cnConfidence} returns a matrix of confidence scores for
	the copy number estimates.

}
\keyword{methods}
\keyword{manip}

