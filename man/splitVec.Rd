\name{splitIndicesByLength}
\alias{splitIndicesByNode}
\alias{splitIndicesByLength}

\title{
  Tools to distribute objects across nodes or by length.
}

\description{
  Tools to distribute objects across nodes or by length.
}

\usage{
splitIndicesByLength(x, lg)
splitIndicesByNode(x)
}

\arguments{
  \item{x}{object to be split}
  \item{lg}{length}
}

\details{
  \code{splitIndicesByLength} splits \code{x} in groups of length \code{lg}.

  \code{splitIndicesByNode} splits \code{x} in N groups (where N is the
  number of compute nodes available).
}

\value{
  List.
}

\author{
  Benilton S Carvalho
}

\seealso{
  split
}
\examples{
x <- 1:100
splitIndicesByLength(x, 8)
splitIndicesByNode(x)
}

\keyword{manip}
