\name{totalCopyNumber}
\alias{totalCopyNumber}
\alias{totalCopyNumber-methods}
\alias{totalCopyNumber,CNSet-method}
\title{
	Compute total copy number 
}
\description{

	Computes total copy number by adding the CA and CB copy number
	estimates at polymorphic loci.  At nonpolymorphic loci, CA is
	returned.

}
\usage{
	totalCopyNumber(object, i, j, ...)
}

\arguments{

  \item{object}{

  Object of class \code{CNSet}.  This object must be annotated with
  information on whether the marker is a SNP or nonpolymorphic.

}

  \item{i}{

  integer vector for subsetting the rows

}

\item{j}{

  integer vector for subsetting the columns

}

\item{...}{

     Ignored.

}

}

\details{

This accessor was created to deal with objects in which the internal
data structure of the \code{assayDataElements} may be either \code{ff}
objects or simple matrices. For \code{ff} or \code{ffdf} objects, we
require specification of \code{i} or \code{j} in order to prevent the
accidental pulling the entire dataset (stored on disk) to memory.

}

\value{
	\code{matrix}
}

\seealso{

	\code{\link{addFeatureAnnotation}}

}

\author{
R. Scharpf
}
\keyword{manip}
\keyword{methods}




