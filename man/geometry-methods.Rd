\name{geometry}
\alias{geometry}
\alias{geometry,DBPDInfo-method}
\alias{geometry,FeatureSet-method}

\title{Array Geometry Information}
\description{
  For a given array, \code{geometry} returns the physical geometry of it.
}
\usage{
geometry(object)
}
\arguments{
  \item{object}{\code{PDInfo} or \code{FeatureSet} object}
}
\examples{
if (require(pd.mapping50k.xba240))
 geometry(pd.mapping50k.xba240)
}
\keyword{manip}
