\name{checkOrder}
\alias{checkOrder}
\alias{chromosomePositionOrder}
\title{Checks whether a eSet-derived class is ordered by chromosome and
  physical position}
\description{
Checks whether a \code{eSet}-derived class (e.g., a \code{SnpSet} or
\code{CNSet} object) is ordered by chromosome and
  physical position
}
\usage{
checkOrder(object, verbose = FALSE)
chromosomePositionOrder(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{SnpSet} or \code{CopyNumberSet}.
  }
  \item{verbose}{
    Logical.
  }
  \item{\dots}{
    additional arguments to \code{order}
  }
}
\details{
  Checks whether the object is ordered by chromosome and physical position.
}
\value{ Logical
}
\author{
R. Scharpf
}


\seealso{
\code{\link{order}}
}
\examples{
data(oligoSetExample)
checkOrder(oligoSet)
oligoSet2 <- chromosomePositionOrder(oligoSet)
checkOrder(oligoSet2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

