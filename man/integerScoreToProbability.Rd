\name{integerScoreToProbability}
\alias{integerScoreToProbability}
\title{
Converts CRLMM integer representation of genotype confidence to probability scale
}
\description{
Converts CRLMM integer representation of genotype confidence to probability scale 
}
\usage{
integerScoreToProbability(x)
}
\arguments{
  \item{x}{integer vector or matrix representation of CRLMM posterior probabilities}
}

\details{

	Posterior probabilities for the CRLMM genotypes, p, are stored
	as integers using the function round(-1000*log2(1-p)).  The
	integer representation of the posterior probabilities can be
	transformed back to the [0,1] scale by 1-2^(-x/1000).

}
\value{
	Numerical value in the interval [0,1].
}
\author{
	R. Scharpf
}
\examples{

}
\keyword{manip}

