\name{i2p}
\alias{i2p}
\alias{p2i}

\title{

	Functions to convert probabilities to integers, or integers to probabilities.

}

\description{

	Probabilities estimated in the \code{crlmm} package are often
	stored as integers to save memory.  We provide a few utility
	functions to go back and forth between the probability and
	integer representations.

}
\usage{
i2p(i)
p2i(p)
}

\arguments{
  \item{i}{  A matrix or vector of integers.}
  \item{p}{  A matrix or vector of probabilities.}
}

\value{

The value returned by \code{i2p} is 

1 - exp(-i/1000)

The value returned by \code{2pi} is 

  as.integer(-1000*log(1-p))

}

\seealso{ \code{\link{confs}}}
\examples{
	i2p(693)
	p2i(0.5)
	i2p(p2i(0.5))
}
\keyword{manip}
