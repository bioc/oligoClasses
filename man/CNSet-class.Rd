\name{CNSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CNSet-class}
\alias{CA,CNSet-method}
\alias{CA<-,CNSet,matrix-method}
\alias{CB,CNSet-method}
\alias{CB<-,CNSet,matrix-method}
\alias{initialize,CNSet-method}

\title{Class "CNSet"}}
\description{A class to store copy number estimates as well as locus-level summaries of the normalized intensities}
\section{Objects from the Class}{

Objects can be created by calls of the form 

\code{new("CNSet", call, CA, CB, callProbability, alleleA, alleleB, phenoData, featureData, experimentData, protocolData, segmentData, emissionPr, position, chromosome, isSnp, annotation, ...)}.

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
Class \code{"\linkS4class{SnpSuperSet}"}, directly.
Class \code{"\linkS4class{AlleleSet}"}, by class "SnpSuperSet", distance 2.
Class \code{"\linkS4class{SnpSet}"}, by class "SnpSuperSet", distance 2.
Class \code{"\linkS4class{eSet}"}, by class "SnpSuperSet", distance 3.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "SnpSuperSet", distance 4.
Class \code{"\linkS4class{Versioned}"}, by class "SnpSuperSet", distance 5.
}
\section{Methods}{
  \describe{
    \item{CA}{\code{signature(object = "CNSet")}: ... }
    \item{CA<-}{\code{signature(object = "CNSet", value = "matrix")}: ... }
    \item{CB}{\code{signature(object = "CNSet")}: ... }
    \item{CB<-}{\code{signature(object = "CNSet", value = "matrix")}: ... }
   \item{initialize}{\code{signature(.Object = "CNSet")}: ... }
   }
}
\author{R. Scharpf}

\seealso{
	\code{\linkS4class{AlleleSet}}, \code{\linkS4class{SnpSuperSet}}
}
\examples{
showClass("CNSet")
}
\keyword{classes}
