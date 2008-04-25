\name{GeneFeatureSet-class}
\docType{class}
\alias{GeneFeatureSet-class}
\title{Class GeneFeatureSet}
\description{A short description}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GeneFeatureSet", assayData, phenoData, featureData, experimentData, annotation, ...)}.
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
  Class \code{\link[Biobase:class.eSet]{eSet}}, by class "FeatureSet",
  distance 2.
  Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, by class "FeatureSet", distance 3.
  Class \code{\link[Biobase:class.Versioned]{Versioned}}, by class "FeatureSet", distance 4.  
}
\section{Methods}{
  No methods defined with class "GeneFeatureSet" in the signature.
}
\author{RS}
\seealso{
}
\examples{
showClass("GeneFeatureSet")
}
\keyword{classes}
