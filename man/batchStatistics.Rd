\name{batchStatistics}
\alias{batchStatistics}
\alias{batchStatistics<-}
\title{

	Accessor for batch statistics uses for copy number estimation
	and storage of model parameters

}

\description{

	The \code{batchStatistics} slot contains statistics estimated
	from each batch that are used to derive copy number estimates.  

}
\usage{
batchStatistics(object)
batchStatistics(object) <- value
}
\arguments{
  \item{object}{  An object of class \code{CNSet}}
  \item{value}{An object of class \code{AssayData}}
}
\details{

	An object of class \code{AssayData} for slot
	\code{batchStatistics} is initialized automatically when
	creating a new \code{CNSet} instance.  Required in the call to
	\code{new} is a factor called \code{batch} whose unique values
	determine the number of columns for each assay data
	element. 

}

\value{

	\code{batchStatics} is an accessor for the slot
	\code{batchStatistics} that returns an object of class
	\code{AssayData}. 

}

\seealso{
	\code{\link{CNSet-class}}, \code{\link{batchNames}}, \code{\link{batch}}
}
\keyword{manip}
\keyword{methods}
