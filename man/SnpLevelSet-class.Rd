\name{SnpLevelSet-class}
\docType{class}
\alias{SnpLevelSet-class}
\alias{calls,SnpLevelSet-method}
\alias{calls<-,SnpLevelSet,matrix-method}
\alias{callsConfidence,SnpLevelSet-method}
\alias{callsConfidence<-,SnpLevelSet,matrix-method}
\alias{chromosome,SnpLevelSet-method}
\alias{chromosome<-,SnpLevelSet,character-method}
\alias{cnConfidence,SnpLevelSet-method}
\alias{cnConfidence<-,SnpLevelSet,matrix-method}
\alias{combine,SnpLevelSet,SnpLevelSet-method}
\alias{copyNumber,SnpLevelSet-method}
\alias{copyNumber<-,SnpLevelSet,matrix-method}
\alias{db,SnpLevelSet-method}
\alias{position,SnpLevelSet-method}

\title{Class "SnpLevelSet"}
\description{Virtual class for SNP-level data}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    
    \item{\code{assayData}:}{The usual restrictions on \code{assayData}:
	contains matrices with equal dimensions, and with column number
	equal to \code{nrow(phenoData)}. Class:\code{\link[Biobase]{class.AssayData}}}
      
    \item{\code{phenoData}:}{See \code{\link[Biobase]{class.eSet}}}
    \item{\code{featureData}:}{See \code{\link[Biobase]{class.eSet}}}
    \item{\code{experimentData}:}{See \code{\link[Biobase]{class.eSet}}}
    \item{\code{annotation}:}{See \code{\link[Biobase]{class.eSet}}}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{eSet}"}, directly.
  Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
  Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}

\section{Methods}{
  \describe{
    \item{calls}{\code{signature(object = "SnpLevelSet")}: Accessor for
  matrix of genotype calls}
    \item{calls<-}{\code{signature(object = "SnpLevelSet", value =
  "matrix")}: Replacement method for genotype calls.}
    \item{callsConfidence}{\code{signature(object = "SnpLevelSet")}:
  Accessor for matrix of confidence scores for the genotype calls}
    \item{callsConfidence<-}{\code{signature(object = "SnpLevelSet",
  value = "matrix")}: Replacement method}
    \item{chromosome}{\code{signature(object = "SnpLevelSet")}:
  See also \code{\link{chromosome}}}
    \item{chromosome<-}{\code{signature(object = "SnpLevelSet", value =
  "character")}: See also \code{\link{chromosome}}}
    \item{cnConfidence}{\code{signature(object = "SnpLevelSet")}:
  Accessor for matrix of confidence scores for the copy number estimates.}
    \item{cnConfidence<-}{\code{signature(object = "SnpLevelSet", value
  = "matrix")}: Replacement method for confidence scores of copy number estimates.}
    \item{combine}{\code{signature(x = "SnpLevelSet", y =
  "SnpLevelSet")}: Useful for combining two objects with the same
  \code{phenoData}. For instance, combine an object created from the
  Hind 50k chip and an object created from the Xba 50k chip.}
    \item{copyNumber}{\code{signature(object = "SnpLevelSet")}: Accessor
  for matrix of copy numbers.}
    \item{copyNumber<-}{\code{signature(object = "SnpLevelSet", value =
  "matrix")}: Replacement method for copy number matrices}
    \item{db}{\code{signature(object = "SnpLevelSet")}: See also \code{\link{db}}}
    \item{position}{\code{signature(object = "SnpLevelSet")}: See also \code{\link{position}}}
	 }
}
\author{RS}
\examples{
showClass("SnpLevelSet")
}
\keyword{classes}
