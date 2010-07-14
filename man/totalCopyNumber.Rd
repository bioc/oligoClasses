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

  Object of class \code{CNSet} 

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

\value{
	\code{matrix}
}

\author{
R. Scharpf
}
\keyword{manip}
\keyword{methods}




