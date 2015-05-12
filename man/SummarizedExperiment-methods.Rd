\name{SummarizedExperiment-methods}
\alias{baf,RangedSummarizedExperiment-method}
\alias{chromosome,RangedSummarizedExperiment-method}
\alias{isSnp,RangedSummarizedExperiment-method}
\alias{lrr,RangedSummarizedExperiment-method}
\alias{baf,SummarizedExperiment-method}
\alias{chromosome,SummarizedExperiment-method}
\alias{isSnp,SummarizedExperiment-method}
\alias{lrr,SummarizedExperiment-method}
\title{Methods for RangedSummarizedExperiment objects}
\description{
  Methods for \link{RangedSummarizedExperiment}.
}
\usage{
\S4method{baf}{RangedSummarizedExperiment}(object)
\S4method{chromosome}{RangedSummarizedExperiment}(object,...)
\S4method{isSnp}{RangedSummarizedExperiment}(object, ...)
\S4method{lrr}{RangedSummarizedExperiment}(object)

}

\arguments{
  \item{object}{
    A \link{RangedSummarizedExperiment} object.
  }
  \item{...}{ignored}
}

\details{

  \code{baf} and \code{lrr} are accessors for the B allele
  frequencies and log R ratio assays (matrices or arrays), respectively,

  \code{chromosome} returns the \code{seqnames} of the \code{rowRanges}.

  \code{isSnp} returns a logical vector for each marker in
  \code{rowRanges} indicating whether the marker targets a SNP
  (nonpolymorphic regions are FALSE).

}

\seealso{
  \code{\linkS4class{RangedSummarizedExperiment}}
}


\keyword{methods}
