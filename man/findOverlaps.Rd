\name{findOverlaps}
\alias{findOverlaps,RangedDataCNV,SnpSet-method}
\alias{findOverlaps,RangedDataCNV,CNSet-method}
\alias{findOverlaps,RangedDataCNV,AnnotatedDataFrame-method}
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

  A matrix. See details.

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

%\examples{
%}
\keyword{methods}
\keyword{manip}
