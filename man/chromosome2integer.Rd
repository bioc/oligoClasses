\name{chromosome2integer}
\alias{chromosome2integer}
\title{Converts chromosome to integer}
\description{
	Coerces character string for chromosome in the pd. annotation packages to integers
}
\usage{
	chromosome2integer(chrom)
}
\arguments{
	\item{chrom}{chromosome}
}
\details{
  This is useful when sorting SNPs in an object by chromosome and
  physical position -- ensures that the sorting is done in the same way
  for different objects.
}
\value{
	  \code{integer}
	  \code{character}
}
\examples{
	chromosome2integer(c(1:22, "X", "Y", "XY", "M"))
}	
\author{R. Scharpf}
\keyword{manip}

