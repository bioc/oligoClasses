\name{SummarizedExperiment-methods}
\alias{baf,SummarizedExperiment-method}
\alias{chromosome,SummarizedExperiment-method}
\alias{isSnp,SummarizedExperiment-method}
\alias{lrr,SummarizedExperiment-method}
\title{Methods for SummarizedExperiment objects}
\description{
  Methods for \link{SummarizedExperiment}.
}
\usage{
\S4method{baf}{SummarizedExperiment}(object)
\S4method{chromosome}{SummarizedExperiment}(object,...)
\S4method{isSnp}{SummarizedExperiment}(object, ...)
\S4method{lrr}{SummarizedExperiment}(object)

}

\arguments{
  \item{object}{
    A \link{SummarizedExperiment} object.
  }
  \item{...}{ignored}
}

\details{

  \code{baf} and \code{lrr} are accessors for the B allele
  frequencies and log R ratio assays (matrices or arrays), respectively,

  \code{chromosome} returns the \code{seqnames} of the \code{rowData}.

  \code{isSnp} returns a logical vector for each marker in
  \code{rowData} indicating whether the marker targets a SNP
  (nonpolymorphic regions are FALSE).

}

\seealso{
  \code{\linkS4class{SummarizedExperiment}}
}


\keyword{methods}
