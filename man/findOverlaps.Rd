\name{findOverlaps}
\alias{findOverlaps,RangedDataCNV,SnpSet-method}
\alias{findOverlaps,RangedDataCNV,CNSet-method}
\alias{findOverlaps,RangedDataCNV,AnnotatedDataFrame-method}
\alias{findOverlaps,RangedDataCNV,RangedDataCNV-method}
\alias{findOverlaps,RangedDataHMM,RangedDataHMM-method}
\alias{findOverlaps}
\title{Find markers that overlap with a query set of ranges}
\description{
  Find indices of markers in a \code{SnpSet} object (the 'subject') that
  overlap with a query set of ranges given by a \code{RangedDataCNV} object.
}
\usage{
findOverlaps(query, subject, maxgap = 0L, minoverlap = 1L, type = c("any", "start", "end", "within", "equal"), select = c("all", "first", "last", "arbitrary"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query}{
    A \code{RangedDataCNV} object.
}
  \item{subject}{
    A \code{SnpSet}, a \code{CNSet} or the \code{featureData} from these
    classes.  The \code{featureData} must be an \code{AnnotatedDataFrame}.
}
  \item{maxgap}{
    Passed to \code{findOverlaps} method for a \code{IRanges} query and
    a \code{IRanges} subject.
}
  \item{minoverlap}{
    Passed to \code{findOverlaps} method for a \code{IRanges} query and
    a \code{IRanges} subject.
}
  \item{type}{
    Passed to \code{findOverlaps} method for a \code{IRanges} query and
    a \code{IRanges} subject.
}
  \item{select}{
    Passed to \code{findOverlaps} method for a \code{IRanges} query and
    a \code{IRanges} subject.
}
  \item{\dots}{
    Passed to \code{findOverlaps} method for a \code{IRanges} query and
    a \code{IRanges} subject.
}
}


\value{

A \code{RangesMatching} object.

}

\details{

    Returns the 'match matrix' from a \code{findOverlaps} query of two
  IRanges objects.  The first column are indices of the ranges in the
  \code{query} (instance of class \code{RangedDataCNV}).  The second
  column are the indices of markers in the \code{subject} (instance of
  class \code{SnpSet}).

}

\author{
R. Scharpf
}


\examples{
	if(require("VanillaICE")){
		data(hmmResults, package="VanillaICE")
		data(oligoSetExample)
		## Find markers in a oligoSnpSet object that overlap with the
		## 2nd range:
		rmatching <- findOverlaps(hmmResults[2, ], oligoSet)
		index.in.range <- matchMatrix(rmatching)[,2]
		features.in.range <- featureNames(oligoSet)[index.in.range]
		##
		## For two RangedDataHMM objects, returns only the ranges
		## that have the same sample name, chromosome, and HMM state
		## A trivial example:
		hmmResults2 <- hmmResults
		hmmResults2$state <- rep(3, nrow(hmmResults))
		findOverlaps(hmmResults, hmmResults2)
	}
}
\keyword{methods}
\keyword{manip}
