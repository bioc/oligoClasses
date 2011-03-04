\name{flags}
\alias{flags}

\title{
	
	Batch-level summary of SNP flags.	

}
\description{
	
	Used to flag SNPs with low minor allele frequencies, or for
	possible problems during the CN estimation step.  Currently,
	this is primarily more for internal use.
	
}
\usage{
flags(object)
}
\arguments{
  \item{object}{
  An object of class \code{CNSet}
}
}

\value{

	A \code{matrix} or \code{ff_matrix} object with rows
	corresponding to markers and columns corresponding to batch.

}

\seealso{
	\code{\link{batchStatistics}}
}
\examples{
x <- matrix(runif(250*96*2, 0, 2), 250, 96*2)
test1 <- new("CNSet", alleleA=x, alleleB=x, call=x, callProbability=x,
	      batch=as.character(rep(letters[1:2], each=96)))
dim(flags(test1))
}
\keyword{manip}
\keyword{methods}
