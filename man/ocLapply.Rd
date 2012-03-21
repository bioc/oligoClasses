\name{ocLapply}
\alias{ocLapply}

\title{
lapply-like function that parallelizes code when possible.
}
\description{
ocLapply is an lapply-like function that checks if ff/snow are loaded
and if the cluster variable is set to execute FUN on a cluster. If these
requirements are not available, then lapply is used.
}
\usage{
ocLapply(X, FUN, ..., neededPkgs)
}

\arguments{
  \item{X}{first argument to FUN.}
  \item{FUN}{function to be executed.}
  \item{\dots}{additional arguments to FUN.}
  \item{neededPkgs}{packages needed to execute FUN on the compute nodes.}
}

\details{
  \code{neededPkgs} is needed when parallel computing is expected to be
  used. These packages are loaded on the compute nodes before the
  execution of FUN.
}

\value{
  A list of length length(X).
}

\author{
  Benilton S Carvalho
}

\seealso{lapply, parStatus}

\keyword{manip}
