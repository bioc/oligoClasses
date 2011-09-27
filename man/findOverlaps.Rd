\name{findOverlaps}
\alias{findOverlaps,RangedDataCNV,SnpSet-method}
\alias{findOverlaps,RangedDataCNV,CNSet-method}
\alias{findOverlaps,RangedDataCNV,AnnotatedDataFrame-method}
\alias{findOverlaps,RangedDataCNV,RangedDataCNV-method}
\alias{findOverlaps,RangedDataHMM,RangedDataHMM-method}
\alias{findOverlaps}
\title{Find markers that overlap with a query set of ranges}
\description{

  Methods for finding overlapping genomic ranges.


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

  When both \code{query} and \code{subject} are
  \code{\linkS4class{RangedDataCNV}} objects, we require that the
  overlapping ranges have the same chromosome and sample id.  If
  \code{query} and \code{subject} are \code{\linkS4class{RangedDataHMM}}
  objects we additionally require that the \code{state} value in the
  \code{\linkS4class{RangedDataHMM}} objects match. For example, if the
  same genomic interval for a subject is called a 'deletion' in
  \code{query} and 'normal' in \code{subject}, the interval is not
  matched. The primary purpose of such queries is to assess the
  concordance of hidden Markov models applied to the same dataset.

  If \code{subject} is a \code{\linkS4class{SnpSet}},
  \code{\linkS4class{CNSet}}, or \code{\linkS4class{AnnotatedDataFrame}}
  with 'chromosome', and 'position' \code{varLabels}, the method returns
  a \code{\linkS4class{RangesMatching}} object indicating which markers
  in \code{subject} overlap with the \code{query} ranges.  This method
  can be useful for finding the set of markers in \code{subject} that
  reside within a given genomic interval in the \code{query}.

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
