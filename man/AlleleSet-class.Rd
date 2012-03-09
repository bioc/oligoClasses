\name{AlleleSet-class}
\Rdversion{1.1}
\docType{class}
\alias{AlleleSet-class}
\alias{AlleleSet}
\alias{allele}
\alias{allele,AlleleSet-method}
\alias{allele,SnpFeatureSet-method}
\alias{bothStrands}
\alias{bothStrands,AlleleSet-method}
\alias{bothStrands,SnpFeatureSet-method}
\alias{db,AlleleSet-method}
\alias{getA,AlleleSet-method}
\alias{getM,AlleleSet-method}

\title{Class "AlleleSet"}
\description{A class for storing the locus-level summaries of the normalized intensities}
\section{Objects from the Class}{
		 Objects can be created by calls of the form \code{new("AlleleSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)}.
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
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{allele}{\code{signature(object = "AlleleSet")}: extract allele
    specific summaries. For 50K (XBA and Hind) and 250K (Sty and
    Nsp) arrays, an additional argument (strand) must be used (allowed
    values: 'sense', 'antisense'.}
    \item{bothStrands}{\code{signature(object = "AlleleSet")}: tests if
      data contains allele summaries on both strands for a given SNP.}
    \item{bothStrands}{\code{signature(object = "SnpFeatureSet")}: tests if
     data contains allele summaries on both strands for a given SnpFeatureSet.}
    \item{db}{\code{signature(object = "AlleleSet")}: link to database connection.}
    \item{getA}{\code{signature(object = "AlleleSet")}: average
      intensities (across alleles)}
    \item{getM}{\code{signature(object = "AlleleSet")}: log-ratio
      (Allele A vs. Allele B)}
  }
}
\author{R. Scharpf}

\seealso{
	\code{\linkS4class{SnpSuperSet}}, \code{\linkS4class{CNSet}}
}
\examples{
showClass("AlleleSet")
## an empty AlleleSet
x <- new("matrix")
new("AlleleSet", senseAlleleA=x, senseAlleleB=x, antisenseAlleleA=x, antisenseAlleleB=x)
##or
new("AlleleSet", alleleA=x, alleleB=x)
}
\keyword{classes}
