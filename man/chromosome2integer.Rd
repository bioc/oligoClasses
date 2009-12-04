\name{chromosome2integer}
\alias{chromosome2integer}
\alias{integer2chromosome}
\title{Converts chromosome to integer or vice versa}
\description{
  Provides a standard way of converting chromosome annotation stored in
  \code{SnpLevelSet} objects to the numeric class.
}
\usage{
chromosome2integer(chrom)
integer2chromosome(chrom)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrom}{chromosome}
}
\details{
  This is useful when sorting SNPs in an object by chromosome and
  physical position -- ensures that the sorting is done in the same way
  for different objects.

  The function currently provides the following conversions:

  "X" -> 23
  "Y" -> 24
  "XY" -> 25  
  "M" -> 26
  }

\value{
  \code{integer}
  \code{character}
}
\author{R. Scharpf}
\keyword{manip}

