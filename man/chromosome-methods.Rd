\name{chromosome-methods}
\docType{methods}
\alias{chromosome-methods}
\alias{chromosome}
\alias{chromosome<-}
\alias{chromosome,AnnotatedDataFrame-method}
\alias{chromosome,SnpSet-method}
\alias{chromosome,gSet-method}
\alias{chromosome,GenomeAnnotatedDataFrame-method}
%\alias{chromosome,RangedDataCNV-method}
\alias{chromosome<-,SnpSet,integer-method}
\alias{chromosome<-,gSet,integer-method}
\alias{chromosome<-,GenomeAnnotatedDataFrame,integer-method}
\alias{chromosome,GRanges-method}
\alias{chromosome,GRangesList-method}
\title{ Methods for function chromosome in package oligoClasses}
\description{
 Methods for function \code{chromosome} in package \pkg{oligoClasses} ~~
}
\section{Methods}{

  The methods for \code{chromosome} extracts the chromosome (represented
  as an integer) for each marker in a \code{eSet}-derived class or a
  \code{AnnotatedDataFrame}-derived class.

  \describe{



\item{\code{signature(object = "AnnotatedDataFrame")}}{
  Accessor for chromosome.
}

\item{\code{signature(object = "eSet")}}{ If 'chromosome' is included in
    \code{fvarLabels(object)}, the integer representation of the
    chromosome will be returned. Otherwise, an error is thrown.  }

  \item{\code{signature(object = "GenomeAnnotatedDataFrame")}}{

  Accessor for chromosome.  If annotation was not available due to a
  missing or non-existent annotation package, the value returned by the
  accessor will be a vector of zero's.}

  \item{}{\code{(chromosome(object) <- value)}: Assign chromosome to the
  \code{AnnotatedDataFrame} slot of an \code{eSet}-derived \code{object}.}

\item{\code{signature(object = "RangedDataCNV")}}{Accessor for chromosome.}
}
}

\note{

  Integer representation: chr X = 23, chr Y = 24, chr XY = 25.  Symbols
  M, Mt, and MT are coded as 26.

}

\examples{
chromosome2integer(c(1:22, "X", "Y", "XY", "M"))
}

\seealso{\code{\link{chromosome2integer}}}
\keyword{methods}
