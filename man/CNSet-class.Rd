\name{CNSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CNSet-class}

\alias{initialize,CNSet-method}
\title{Class "CNSet"}
\description{A class to store parameters for CN estimation as well as locus-level summaries of the normalized intensities}
\section{Objects from the Class}{

		 An object from the class is not generally intended to
		 be initialized by the user, but returned by the
		 \code{genotype} function.  The example below
		 initializes an empty CNSet object.

}
\section{Slots}{
  \describe{
	  \item{\code{lM}:}{Object of class \code{"list_or_ffdf"}. For storing parameters related to copy number estimation.}
	  \item{\code{assayData}:}{Object of class \code{"AssayData"} }
	  \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} }
	  \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} }
	  \item{\code{experimentData}:}{Object of class \code{"MIAME"} }
	  \item{\code{annotation}:}{Object of class \code{"character"}}
	  \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} }
	  \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SnpSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "SnpSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "SnpSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "SnpSet", distance 4.
}
\section{Methods}{
  \describe{
  \item{A}{\code{signature(object = "CNSet")}: accessor for normalized intensities for A allele}
  \item{B}{\code{signature(object = "CNSet")}: accessor for normalized intensities for B allele }
  \item{calls}{\code{signature(object = "CNSet")}: accessor for crlmm genotype calls}
  \item{confs}{\code{signature(object = "CNSet")}: accessor for genotype confidence scores}
}
}
\author{R. Scharpf}

\seealso{
	\code{\linkS4class{AlleleSet}}, \code{\linkS4class{SnpSuperSet}}
}
\examples{
	##An empty CNSet
	x <- new("matrix")
	new("CNSet", alleleA=x, alleleB=x, call=x, callProbability=x, lM=list())
}
	
\keyword{classes}
