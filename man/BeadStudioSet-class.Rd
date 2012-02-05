\name{BeadStudioSet-class}
\Rdversion{1.1}
\docType{class}
\alias{BeadStudioSet-class}
\alias{baf,BeadStudioSet-method}
\alias{baf<-,BeadStudioSet-method}
\alias{baf<-}
\alias{copyNumber,BeadStudioSet-method}
\alias{copyNumber<-,BeadStudioSet,ANY-method}
\alias{coerce,BeadStudioSet,data.frame-method}
\alias{initialize,BeadStudioSet-method}
\alias{lrr,BeadStudioSet-method}
\alias{lrr<-,BeadStudioSet-method}
\alias{lrr<-}

\title{Class \code{"BeadStudioSet"}}
\description{
  A container for log R ratios and B allele frequencies from SNP arrays.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BeadStudioSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, baf, lrr, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{featureData}:}{Object of class \code{"GenomeAnnotatedDataFrame"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{gSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "gSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "gSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "gSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{baf}{\code{signature(object = "BeadStudioSet")}: accessor for
      B allele frequencies}
    \item{}{\code{baf(object) <- value} replacement
    method for B allele frequencies: \code{value} must be a matrix.}
    \item{coerce}{\code{signature(from = "BeadStudioSet", to =
	"data.frame")}: coerce to data.frame with column headers 'lrr',
      'baf', 'x' (physical position with unit Mb), 'id', and 'is.snp'.
      Used for plotting with lattice.}
    \item{}{\code{copyNumber(object)}:  accessor for the log R ratios.}
    \item{}{\code{copyNumber(object) <- value}:  replacement method for
      the log R ratios}
    \item{initialize}{\code{signature(.Object = "BeadStudioSet")}:
      constructs an instance of the class}
    \item{lrr}{\code{signature(object = "BeadStudioSet")}: accessor for
      log R ratios}
    \item{}{\code{lrr(object) <- value} replacement
    method for log R ratios: \code{value} must be a matrix.}
}
}

\author{
R. Scharpf
}

\examples{
new("BeadStudioSet")
}
\keyword{classes}
