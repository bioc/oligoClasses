\name{batch}
\alias{batch}
\alias{batchNames}
\alias{batchNames<-}

\title{
	The batch variable for the samples.
}

\description{

Copy number estimates are susceptible to systematic differences between
  groups of samples that were processed at different times or by
  different labs.  While 'batch' is often unknown, a useful surrogates
  is often the scan date of the arrays (e.g., the month of the calendar
  year) or the 96 well chemistry plate on which the samples were arrayed
  during lab processing.

}

\usage{
batch(object)
batchNames(object)
batchNames(object) <- value
}

\arguments{
  \item{object}{
  An object of class \code{CNSet}.
}

  \item{value}{

  For 'batchNames', the value must be a character string corresponding
  of the unique batch names.

}

}

\value{

  The method 'batch' returns a \code{character} vector that has the same
  length as the number of samples in the \code{CNSet} object.

}

\author{
R. Scharpf
}

\seealso{
\code{\link{CNSet-class}}
}
\examples{
a <- matrix(1:25, 5, 5)
colnames(a) <- letters[1:5]
object <- new("CNSet", alleleA=a, batch=rep("batch1", 5))
batch(object)
batchNames(object)
}
\keyword{methods}

