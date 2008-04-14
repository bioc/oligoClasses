\name{chromosome}
\alias{chromosome}
\alias{chromosome<-}
\alias{chromosome<-,AnnotatedDataFrame,character-method}
\title{Accessor for chromosome information}
\description{
  Returns chromosome information.
}
\usage{
chromosome(object)
chromosome(object) <- value
}
\arguments{
  \item{object}{object inheriting from \code{SnpLevelSet}}
  \item{value}{vector of class \code{character}}
}
\value{
  Vector with chromosome information.  Chromosomes are recorded as
  character strings (e.g., "1", "2", "X", ...)
}
\author{R. Scharpf}
\keyword{manip}
