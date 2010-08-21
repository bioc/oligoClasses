\name{batch}
\alias{batch}
\alias{batchNames}
\alias{batchNames<-}

\title{

	The batch variable for the samples.

}

\description{

Copy number estimates are susceptible to systematic differences
  between groups of samples that were processed at different times or
  by different labs.  Analysis algorithms that do not adjust for batch
  effects are prone to spurious measures of association.  While
  'batch' is often unknown, a useful surrogates are the scan date of
  the arrays or the 96 well chemistry plate on which the samples were
  arrayed during lab processing.

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

\details{

}

\value{

	The method 'batch' returns a factor that has the same length
	as the number of samples in the \code{CNSet} object.

	The method 'batchNames' returns the unique batches as a
	character string.  The batch labels for each element in the
	\code{LinearModelParameter} class can be reassigned using the
	'batchNames<-' replacement method.

}

\author{
R. Scharpf

}

\seealso{
\code{\link{CNSet-class}}
}
\examples{

x <- matrix(runif(250*96*2, 0, 2), 250, 96*2)
test1 <- new("CNSet", alleleA=x, alleleB=x, call=x, callProbability=x,
	     batch=as.factor(rep(letters[1:2], each=96)))
batchNames(test1) ##unique batches
batch(test1)
test1[1:20, 1:10]
##just NA's
nu(test1, "A")

## using ff objects
library(ff)
x2 <- initializeBigMatrix("smallx", nr=250, nc=96*2)
x2[,] <- as.numeric(x)
test2 <- new("CNSet", alleleA=x, alleleB=x, call=x, callProbability=x, batch=as.factor(rep(letters[1:2], each=96)))
test2
batchNames(test2) ##unique batches
batch(test2)
(test2.sub <- test2[1:20, 1:10])
##just NA's
nu(test2.sub, "A")

}
\keyword{methods}

