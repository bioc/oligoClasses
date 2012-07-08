\name{chromosome2integer}
\alias{chromosome2integer}
\alias{integer2chromosome}
\title{Converts chromosome to integer}
\description{
	Coerces character string for chromosome in the pd. annotation packages to integers
}
\usage{
	chromosome2integer(chrom)
        integer2chromosome(intChrom)
}
\arguments{
  \item{chrom}{A one or 2 letter character string (e.g, "1", "X", "Y",
    "MT", "XY") }
  \item{intChrom}{An integer vector with values 1-25 possible}
}
\details{
  This is useful when sorting SNPs in an object by chromosome and
  physical position -- ensures that the sorting is done in the same way
  for different objects.
}
\value{
  \code{integer2chromosome} returns a vector of character string
  indicating the chromosome the same length
  as \code{intChrom}.  \code{chromosome2integer} returns a vector of
  integers the same length as the number of elements in the \code{chrom} vector.
}
\examples{
	chromosome2integer(c(1:22, "X", "Y", "XY", "M"))
	integer2chromosome(chromosome2integer(c(1:22, "X", "Y", "XY", "M")))
}
\author{R. Scharpf}
\keyword{manip}

