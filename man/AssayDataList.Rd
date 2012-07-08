\name{AssayDataList}
\alias{AssayDataList}
\title{Create a list of assay data elements}
\description{The eSetList-derived classes have an assayDataList slot
  instead of an assayData slot.
}
\usage{
AssayDataList(storage.mode = c("lockedEnvironment", "environment", "list"), ...)
}

\arguments{
  \item{storage.mode}{See \code{assayDataNew}.}
  \item{\dots}{Named lists of matrices}
}

\value{
  environment
}
\author{
R.Scharpf
}

\seealso{
  \code{\link{assayDataNew}}
}
\examples{
r <- replicate(5, matrix(rnorm(25),5,5), simplify=FALSE)
r <- lapply(r, function(x,dns) {dimnames(x) <- dns; return(x)}, dns=list(letters[1:5], LETTERS[1:5]))
ad <- AssayDataList(r=r)
ls(ad)
}
\keyword{manip}

