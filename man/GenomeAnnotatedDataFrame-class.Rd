\name{GenomeAnnotatedDataFrame-class}
\Rdversion{1.1}
\docType{class}
\alias{GenomeAnnotatedDataFrame-class}
%\alias{chromosome,GenomeAnnotatedDataFrame-method}
%\alias{chromosome<-,GenomeAnnotatedDataFrame,integer-method}
\alias{coerce,AnnotatedDataFrame,GenomeAnnotatedDataFrame-method}
\alias{initialize,GenomeAnnotatedDataFrame-method}
%\alias{isSnp,GenomeAnnotatedDataFrame,ANY-method}
%\alias{position,GenomeAnnotatedDataFrame-method}
\alias{updateObject,AnnotatedDataFrame-method}

\title{Class \code{"GenomeAnnotatedDataFrame"}}
\description{
  AnnotatedDataFrame with genomic coordinates (chromosome, position)
}

\section{Slots}{
  \describe{
    \item{\code{varMetadata}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{dimLabels}:}{Object of class \code{"character"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}

\section{Extends}{
Class \code{"\linkS4class{AnnotatedDataFrame}"}, directly.
Class \code{"\linkS4class{Versioned}"}, by class "AnnotatedDataFrame", distance 2.
}

\section{Methods}{
  \describe{
    \item{chromosome}{\code{signature(object = "GenomeAnnotatedDataFrame")}: the chromosome (coded as integers) }
    \item{chromosome<-}{\code{signature(object = "GenomeAnnotatedDataFrame", value = "integer")}: replacement method for chromosome}
    \item{coerce}{\code{signature(from = "AnnotatedDataFrame", to = "GenomeAnnotatedDataFrame")}: coerce an AnnotatedDataFrame with chromosome, position, and isSnp varLabels to a GenomeAnnotatedDataFrame }
    \item{initialize}{\code{signature(.Object = "GenomeAnnotatedDataFrame")}: ... }
    \item{isSnp}{\code{signature(object = "GenomeAnnotatedDataFrame", pkgname = "ANY")}: Many platforms include polymorphic and nonpolymorphic markers.  isSnp evalutes to \code{TRUE} if the marker is polymorphic.}
    \item{position}{\code{signature(object = "GenomeAnnotatedDataFrame")}: Physical position in the genome}
}
}

\author{
R. Scharpf
}

\seealso{
  \code{\linkS4class{AnnotatedDataFrame}}
}

\examples{
new("GenomeAnnotatedDataFrame")
if(require("pd.mapping50k.hind240") && require("pd.mapping50k.xba240")){
	data(locusLevelData)
	object <- GenomeAnnotatedDataFrameFrom(locusLevelData[["genotypes"]],
					       annotationPkg=locusLevelData[["platform"]])
}
}
\keyword{classes}
