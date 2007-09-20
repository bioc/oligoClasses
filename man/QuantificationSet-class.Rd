\name{QuantificationSet-class}
\docType{class}
\alias{QuantificationSet-class}

\title{Class "QuantificationSet"}
\description{A virtual class to store summarized measures.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
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
Class \code{\link[Biobase:class.eSet]{eSet}}, directly.
Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, by class "eSet", distance 2.
Class \code{\link[Biobase:class.Versioned]{Versioned}}, by class "eSet", distance 3.
}
\section{Methods}{
No methods defined with class "QuantificationSet" in the signature.
}
\seealso{
  \code{\linkS4class{SnpQSet}}, \code{\linkS4class{SnpCnvQSet}}
}
\examples{
showClass("QuantificationSet")
}
\keyword{classes}
