\name{SnpCnvFeatureSet-class}
\docType{class}
\alias{SnpCnvFeatureSet-class}

\title{Class "SnpCnvFeatureSet"}
\description{A class to store feature-level intensities from arrays that
contain both SNP and CNV probes.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SnpCnvFeatureSet", assayData, phenoData, featureData, experimentData, annotation, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{manufacturer}:}{Object of class \code{"character"}}
    \item{\code{platform}:}{Object of class \code{"character"}}
    \item{\code{assayData}:}{Object of class \code{"AssayData"}}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}}
    \item{\code{annotation}:}{Object of class \code{"character"}}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SnpFeatureSet}"}, directly.
Class \code{"\linkS4class{FeatureSet}"}, by class "SnpFeatureSet", distance 2.
Class \code{"\linkS4class{eSet}"}, by class "SnpFeatureSet", distance 3.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "SnpFeatureSet", distance 4.
Class \code{"\linkS4class{Versioned}"}, by class "SnpFeatureSet", distance 5.
}
\section{Methods}{
No methods defined with class "SnpCnvFeatureSet" in the signature.
}
\examples{
showClass("SnpCnvFeatureSet")
}
\keyword{classes}
