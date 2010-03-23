\name{requireClusterPkgSet}
\alias{requireClusterPkg}
\alias{requireClusterPkgSet}

\title{
  Package loaders for clusters.
}
\description{
  Package loaders for clusters.
}

\usage{
requireClusterPkgSet(packages)
requireClusterPkg(...)
}

\arguments{
  \item{packages}{character vector with the names of the packages to be
    loaded on the compute nodes.}
  \item{\dots}{same arguments as \code{require}}
}
\details{
  \code{requireClusterPkgSet} applies \code{require} for a set of
  packages on the cluster nodes.

  \code{requireClusterPkg} applies \code{require} for *ONE* package on
  the cluster nodes and accepts every argument taken by \code{require}.
}
\value{
  Logical.
}

\author{
  Benilton S Carvalho
}

\seealso{require}

\keyword{manip}
