\name{featuresInRange}
\alias{featuresInRange}
\title{
  Find feature index in a genomic interval
}
\description{
    Find feature index of a \code{SnpSet} object in a
    \code{RangedDataCNV} object.
}
\usage{
featuresInRange(object, range, FRAME = 0, FRAME.LEFT, FRAME.RIGHT, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{SnpSet}.
}
  \item{range}{
    A \code{RangedDataCNV} object.
}
\item{FRAME}{
  Integer (basepairs).  The distance from start and end coordinates of
  the genomic interval. Useful for retrieving indices that 'frame' the
  genomic interval of interest.
}
  \item{FRAME.LEFT}{
    Integer.  If missing, set equal to FRAME.
}
  \item{FRAME.RIGHT}{
    Integer.  If missing, set equal to FRAME.
}
  \item{\dots}{
    Ignored
}
}
\details{
  Extract marker indices in \code{object} that occur within a genomic interval.
}
\value{
  Vector of integers.
}
\author{
R. Scharpf
}

%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

