\name{requireClusterPkgSet}
\alias{requireClusterPkg}
\alias{requireClusterPkgSet}
\alias{requireClusterPkg-deprecated}
\alias{requireClusterPkgSet-deprecated}

\title{
  DEPRECATED FUNCTIONS.
  Package loaders for clusters.
}
\description{
  Package loaders for clusters.
}

\usage{
requireClusterPkgSet(packages)
requireClusterPkg(pkg, character.only)
}

\arguments{
  \item{packages}{character vector with the names of the packages to be
    loaded on the compute nodes.}
  \item{pkg}{name of a package given as a name or literal character string}	
  \item{character.only}{a logical indicating whether `pkg' can be
          assumed to be a character string}
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
