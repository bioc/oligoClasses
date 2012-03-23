\name{integerMatrix}
\alias{integerMatrix}
\alias{integerArray}
\title{
  Coerce numeric matrix (or array) to a matrix (array) of integers, retaining dimnames.
}
\description{
  Coerce numeric matrix to matrix of integers, retaining dimnames.
}
\usage{
integerMatrix(x, scale = 100)
integerArray(x, scale=100)
}
\arguments{
  \item{x}{a \code{matrix} or \code{array}}
  \item{scale}{scalar (numeric).  If not 1, \code{x} is multiplied by
    \code{scale} prior to coercing to a matrix of integers.}
}

\value{
  A \code{matrix} or \code{array} of integers.
}

\author{
R. Scharpf
}


\examples{
x <- matrix(rnorm(10), 5, 2)
rownames(x) = letters[1:5]
i <- integerMatrix(x, scale=100)
}
\keyword{manip}
