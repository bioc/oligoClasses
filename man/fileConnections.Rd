\name{fileConnections}
\alias{open}
\alias{openff}
\alias{openff,CNSet-method}
\alias{close}
\alias{closeff}
\alias{closeff,CNSet-method}
\alias{open,matrix-method}
\alias{open,numeric-method}
\alias{open,array-method}
\alias{close,matrix-method}
\alias{close,numeric-method}
\alias{close,array-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Open and close methods for matrices and numeric vectors
}
\description{

  \code{CNSet} objects can contain \code{ff}-derived objects that
  contain pointers to files on disk, or ordinary matrices.  Here we
  define open and close methods for ordinary matrices and vectors that
  that simply pass back the original matrix/vector.

}
\usage{
open(con, ...)
openff(object)
closeff(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{
    matrix or vector
  }
  \item{object}{ A \code{CNSet} object.}
  \item{\dots}{Ignored}
}

\value{
  not applicable
}
\author{
R. Scharpf
}

\examples{
open(rnorm(15))
open(matrix(rnorm(15), 5,3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
