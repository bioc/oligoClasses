\name{SnpSuperSet-class}
\Rdversion{1.1}
\docType{class}
\alias{SnpSuperSet-class}
\alias{SnpSuperSet}
\alias{initialize,SnpSuperSet-method}

\title{Class "SnpSuperSet"}
\description{A class to store locus-level summaries of the quantile normalized intensities, genotype calls, and genotype confidence scores}
\section{Objects from the Class}{
\code{new("SnpSuperSet", allelea=alleleA, alleleB=alleleB, call=call, callProbability, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{AlleleSet}"}, directly.
Class \code{"\linkS4class{SnpSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "AlleleSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "AlleleSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "AlleleSet", distance 4.
}
\section{Methods}{
No methods defined with class "SnpSuperSet" in the signature.
}
\author{R. Scharpf}
\seealso{
	\code{\linkS4class{AlleleSet}}
}
\examples{
showClass("SnpSuperSet")
## empty object from the class
x <- new("matrix")
new("SnpSuperSet", alleleA=x, alleleB=x, call=x, callProbability=x)
}
\keyword{classes}
