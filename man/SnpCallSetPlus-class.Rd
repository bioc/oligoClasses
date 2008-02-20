\name{SnpCallSetPlus-class}
\docType{class}
\alias{SnpCallSetPlus-class}
\alias{SnpCnvCallSetPlus-class}
\alias{calculateCopyNumber,SnpCallSetPlus-method}
\alias{initialize,SnpCallSetPlus-method}
\alias{initialize,SnpCnvCallSetPlus-method}
\alias{coerce,SnpCallSetPlus,oligoSnpSet-method}

\title{Class "SnpCallSetPlus"}
\description{Description for Genotyping}
\section{Objects from the Class}{
  Objects of this class are unions of Quantification and Calls.
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
Class \code{"\linkS4class{SnpQSet}"}, directly.
Class \code{"\linkS4class{SnpCnvQSet}"}, directly.
Class \code{"\linkS4class{SnpCallSet}"}, directly.
Class \code{"\linkS4class{QuantificationSet}"}, by class "SnpQSet", distance 2.
Class \code{\link[Biobase:class.eSet]{eSet}}, by class "SnpQSet", distance 3.
Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, by class "SnpQSet", distance 4.
Class \code{\link[Biobase:class.Versioned]{Versioned}}, by class "SnpQSet", distance 5.
Class \code{"\linkS4class{QuantificationSet}"}, by class "SnpCnvQSet", distance 2.
Class \code{\link[Biobase:class.eSet]{eSet}}, by class "SnpCnvQSet", distance 3.
Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, by class "SnpCnvQSet", distance 4.
Class \code{\link[Biobase:class.Versioned]{Versioned}}, by class "SnpCnvQSet", distance 5.
Class \code{"\linkS4class{SnpLevelSet}"}, by class "SnpCallSet", distance 2.
Class \code{\link[Biobase:class.eSet]{eSet}}, by class "SnpCallSet", distance 3.
Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, by class "SnpCallSet", distance 4.
Class \code{\link[Biobase:class.Versioned]{Versioned}}, by class "SnpCallSet", distance 5.
}
\section{Methods}{
  \describe{
    \item{calculateCopyNumber}{\code{signature(object="SnpCallSetPlus")}
      See \code{\link{calculateCopyNumber}}}
    }
}

\examples{
showClass("SnpCallSetPlus")
}
\keyword{classes}
