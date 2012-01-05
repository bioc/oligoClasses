\name{gSet-class}
\Rdversion{1.1}
\docType{class}
\alias{gSet-class}
\alias{checkOrder,gSet-method}
\alias{db,gSet-method}

\title{Container for objects with genomic annotation on SNPs}
\description{
  Container for objects with genomic annotation on SNPs
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{featureData}:}{Object of class \code{"GenomeAnnotatedDataFrame"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
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
    \item{checkOrder}{\code{signature(object = "gSet")}: ... }
    \item{chromosome}{\code{signature(object = "gSet")}: ... }
    \item{chromosome<-}{\code{signature(object = "gSet", value = "integer")}: ... }
    \item{db}{\code{signature(object = "gSet")}: ... }
    \item{isSnp}{\code{signature(object = "gSet")}: ... }
    \item{position}{\code{signature(object = "gSet")}: ... }
	 }
}

\author{
R. Scharpf
}


\seealso{
\code{\link{chromosome}}, \code{\link{position}}, \code{\link{isSnp}}
}
\examples{
showClass("gSet")
}
\keyword{classes}
