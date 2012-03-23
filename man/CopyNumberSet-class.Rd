\name{CopyNumberSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CopyNumberSet-class}
\alias{CopyNumberSet}
\alias{cnConfidence,CopyNumberSet-method}
\alias{cnConfidence<-,CopyNumberSet,matrix-method}
\alias{copyNumber,CopyNumberSet-method}
\alias{copyNumber<-,CopyNumberSet,matrix-method}
\alias{initialize,CopyNumberSet-method}
\alias{checkOrder,CopyNumberSet-method}

\title{Class \code{"CopyNumberSet"}}
\description{
  Container for storing total copy number estimates and confidence
  scores of the copy number estimates.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CopyNumberSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, copyNumber, cnConfidence, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{cnConfidence}{\code{signature(object = "CopyNumberSet")}: ... }
    \item{cnConfidence<-}{\code{signature(object = "CopyNumberSet", value = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "CNSet", to = "CopyNumberSet")}: ... }
    \item{copyNumber}{\code{signature(object = "CopyNumberSet")}: ... }
    \item{copyNumber<-}{\code{signature(object = "CopyNumberSet", value = "matrix")}: ... }
    \item{initialize}{\code{signature(.Object = "CopyNumberSet")}: ... }
	 }
}
\author{
  R. Scharpf
}

\note{

  This container is primarily for platforms for which genotypes are
  unavailable.  As \code{oligoSnpSet} extends this class, methods
  related to total copy number that do not depend on genotypes can be
  defined at this level.


}

\seealso{

  For genotyping platforms, total copy number estimates and genotype
  calls can be stored in the \code{\linkS4class{oligoSnpSet}} class.

}
\examples{
showClass("CopyNumberSet")
cnset <- new("CopyNumberSet")
ls(Biobase::assayData(cnset))
}
\keyword{classes}
