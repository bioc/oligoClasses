\name{SnpCnvQSet-class}
\docType{class}
\alias{SnpCnvQSet-class}

\title{Class "SnpCnvQSet"}
\description{A class to store summarized measures from
  \code{SnpCnvFeatureSet} objects}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SnpCnvQSet", assayData, phenoData, featureData, experimentData, annotation, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"}}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}}
    \item{\code{annotation}:}{Object of class \code{"character"}}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
  }
}
\section{Extends}{
Class \code{\linkS4class{QuantificationSet}}, directly.
Class \code{\link[Biobase:class.eSet]{eSet}}, by class "QuantificationSet", distance 2.
Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, by class "QuantificationSet", distance 3.
Class \code{\link[Biobase:class.Versioned]{Versioned}}, by class "QuantificationSet", distance 4.
}
\section{Methods}{
No methods defined with class "SnpCnvQSet" in the signature.
}
\examples{
showClass("SnpCnvQSet")
}
\keyword{classes}
