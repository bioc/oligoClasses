\name{SnpSet-methods}
%\docType{methods}
\alias{calls}
\alias{calls,SnpSet-method}
\alias{calls<-}
\alias{calls<-,SnpSet,matrix-method}
\alias{confs}
\alias{confs,SnpSet-method}
\alias{confs<-}
\alias{confs<-,SnpSet,matrix-method}
\alias{checkOrder,SnpSet-method}

\title{Accessors and methods for SnpSet objects}

\description{
  Utility functions for accessing data in \code{SnpSet} objects.
}

\usage{
calls(object)
calls(object) <- value
confs(object, transform=TRUE)
confs(object) <- value
}

\arguments{
  \item{object}{A \link{SnpSet} object.}
  \item{transform}{Logical. Whether to transform the integer
  representation of the confidence score (for memory efficiency) to a
  probability. See details.}
\item{value}{A matrix.}
}

\details{

  \code{calls} returns the genotype calls. CRLMM stores genotype calls
  as integers (1 - AA; 2 - AB; 3 - BB).

  \code{confs} returns the confidences associated with the genotype
  calls. The current implementation of CRLMM stores the confidences as
  integers to save memory on disk by using the transformation:

  round(-1000*log2(1-p)),

  where 'p' is the posterior probability of the call.  \code{confs} is
  a convenience function that transforms the integer representation
  back to a probability.  Note that if the assayData elements of the
  \code{SnpSet} objects are \code{ff_matrix} or \code{ffdf}, the
  \code{confs} function will return a warning. For such objects, one
  should first subset the \code{ff} object and coerce to a matrix,
  then apply the above conversion.  The function
  \code{snpCallProbability} for the \code{callProbability} slot of
  \code{SnpSet} objects.  See the examples below.

  \code{checkOrder} checks whether the object is ordered by chromosome
  and physical position, evaluating to TRUE or FALSE.

}

\note{
  Note that the replacement method for \code{confs<-} expects a matrix
  of probabilities and will automatically convert the probabilities to
  an integer representation.  See details for the conversion.

  The accessor \code{\link{snpCallProbability}} is an accessor for the
  'callProbability' element of the \code{assayData}.  The name can be
  misleading, however, as the accessor will not return a probability if
  the call probabilities are represented as integers.

}

%\section{Methods}{
%	\describe{
%	\item{\code{initialize(SnpSet)}:}{Object instantiation, used by
%	  \code{new}; not to be called directly by the user.}
%	\item{\code{calls(object)}:}{accessor for genotype calls}
%	\item{\code{confs(object)}:}{accessor for crlmm genotype confidence scores}
%	}
%}
\seealso{

  The helper functions \code{\link{p2i}} converts probabilities to
  integers and \code{\link{i2p}} converts integers to probabilities.

  See \code{\link{order}} and \code{\link{checkOrder}}.

}

\examples{

	theCalls <- matrix(sample(1:3, 20, rep=TRUE), nc=2)
	p <- matrix(runif(20), nc=2)
	integerRepresentation <- matrix(as.integer(round(-1000*log(1-p))), 10, 2)
	obj <- new("SnpSet2", call=theCalls, callProbability=integerRepresentation)
	calls(obj)
        confs(obj)  ## coerces to probability scale
        int <- Biobase::snpCallProbability(obj) ## not necessarily a probability
        p3 <- i2p(int)  ## to convert back to a probability

}
\keyword{manip}
