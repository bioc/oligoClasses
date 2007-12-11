\name{position}
\alias{position}
\alias{position,SnpLevelSet-method}
\title{Accessor to position information}
\description{
  \code{position} will return the genomic position of a SNP.
}
\usage{
position(object)
}
\arguments{
  \item{object}{object inheriting from \code{SnpLevelSet}}
}
\details{
  \code{position} will return genomic position of a SNP (number of
  basepairs from the 5-prime chromosomal end)
}
\value{
  an integer
}
\author{R. Scharpf}
\keyword{manip}
