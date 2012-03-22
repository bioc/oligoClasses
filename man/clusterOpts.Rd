\name{ocSamples}
\alias{ocSamples}
\alias{ocProbesets}

\title{
  Cluster and large dataset management utilities.
}
\description{
  Tools to simplify management of clusters via 'snow' package and large
  dataset handling through the 'bigmemory' package.
}
\usage{
ocSamples(n)
ocProbesets(n)
}

\arguments{
  \item{n}{integer representing the maximum number of samples/probesets
    to be processed simultaneously on a compute node.}
}

\details{
  Some methods in the oligo/crlmm packages, like \code{backgroundCorrect},
  \code{normalize}, \code{summarize} and \code{rma} can use a cluster
  (set through the 'foreach' package). The use of cluster features is
  conditioned on the availability of the 'ff' (used to
  provide shared objects across compute nodes) and 'foreach' packages.

  To use a cluster, 'oligo/crlmm' checks for three requirements: 1) 'ff'
  is loaded; 2) an adaptor for the parallel backend (like 'doMPI',
  'doSNOW', 'doMC') is loaded and registered.

  If only the 'ff' package is available and loaded (in addition
  to the caller package - 'oligo' or 'crlmm'), these methods will allow
  the user to analyze datasets that would not fit in RAM at the expense
  of performance.

  In the situations above (large datasets and cluster), oligo/crlmm uses the
  options \code{ocSamples} and \code{ocProbesets} to limit the
  amount of RAM used by the machine(s). For example, if ocSamples is
  set to 100, steps like background correction and normalization process
  (in RAM) 100 samples simultaneously on each compute node. If
  ocProbesets is set to 10K, then summarization processes 10K
  probesets at a time on each machine.
}

\section{Warning}{
  In both scenarios (large dataset and/or cluster use), there
  is a penalty in performance because data are written to disk (to
  either minimize memory footprint or share data across compute nodes).
}

\author{
  Benilton Carvalho
}

\examples{
   if(require(doMC)) {
       registerDoMC()
       ## tasks like summarize()
   }
}
\keyword{manip}

