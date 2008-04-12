\name{GeneSTFeatureSet-class}
\docType{class}
\alias{GeneSTFeatureSet-class}
\title{Class GeneSTFeatureSet}
\description{A short description}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GeneSTFeatureSet", assayData, phenoData, featureData, experimentData, annotation, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{manufacturer}:}{Object of class \code{"character"}}
    \item{\code{platform}:}{Object of class \code{"character"} }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{FeatureSet}"}, directly.
  Class \code{"\linkS4class{eSet}"}, by class "FeatureSet", distance 2.
  Class \code{"\linkS4class{VersionedBiobase}"}, by class "FeatureSet", distance 3.
  Class \code{"\linkS4class{Versioned}"}, by class "FeatureSet", distance 4.
}
\section{Methods}{
  No methods defined with class "GeneSTFeatureSet" in the signature.
}
\author{RS}
\seealso{
}
\examples{
showClass("GeneSTFeatureSet")
}
\keyword{classes}
