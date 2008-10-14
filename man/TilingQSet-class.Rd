\name{TilingQSet-class}
\docType{class}
\alias{TilingQSet-class}
\alias{initialize,TilingQSet-method}

\title{Class "TilingQSet"}
\description{A place-holder for future uses of Quantification on Tiling Arrays.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TilingQSet", assayData, M, phenoData, featureData, experimentalData, annotation)}.
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
Class \code{"\linkS4class{QuantificationSet}"}, directly.
Class \code{"\linkS4class[Biobase]{eSet}"}, by class "QuantificationSet", distance 2.
Class \code{"\linkS4class[Biobase]{VersionedBiobase}"}, by class "QuantificationSet", distance 3.
Class \code{"\linkS4class[Biobase]{Versioned}"}, by class "QuantificationSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "TilingQSet")}: ... }
	 }
}
\author{Benilton Carvalho}
\note{Alpha-version}
\examples{
showClass("TilingQSet")
}
\keyword{classes}
