\name{ldSetOptions}
\alias{ldSetOptions}
\alias{ldStatus}
\alias{ldPath}

\title{
  Set/check large dataset options.
}
\description{
  Set/check large dataset options.
}
\usage{
ldSetOptions(nsamples=100, nprobesets=20000, path=getwd(), verbose=FALSE)
ldStatus(verbose=FALSE)
ldPath(path)
}

\arguments{
  \item{nsamples}{number of samples to be processed at once.}
  \item{nprobesets}{number of probesets to be processed at once.}
  \item{path}{path where to store large dataset objects.}
  \item{verbose}{verbosity (logical).}
}

\details{
  Some functions in oligo/crlmm can process data in batches to minimize
  memory footprint. When using this feature, the 'ff' package resources
  are used (and possibly combined with cluster resources set in
  options() via 'snow' package).

  Methods that are executed on a sample-by-sample manner can use
  ocSamples() to automatically define how many samples are processed at
  once (on a compute node). Similarly, methods applied to probesets can
  use ocProbesets(). Users should set these options appropriately.

  \code{ldStatus} checks the support for large datasets.

  \code{ldPath} checks where ff files are stored.
}

\author{
  Benilton S Carvalho
}

\seealso{ocSamples, ocProbesets}

\examples{
ldStatus(TRUE)
}

\keyword{manip}
