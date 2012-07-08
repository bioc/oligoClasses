\name{GenomeAnnotatedDataFrame-class}
\Rdversion{1.1}
\docType{class}
\alias{GenomeAnnotatedDataFrame-class}
\alias{GenomeAnnotatedDataFrame}
\alias{coerce,AnnotatedDataFrame,GenomeAnnotatedDataFrame-method}
\alias{initialize,GenomeAnnotatedDataFrame-method}
\alias{makeFeatureRanges,GenomeAnnotatedDataFrame-method}
\alias{position<-}
\alias{position<-,GenomeAnnotatedDataFrame,integer-method}
\alias{updateObject,GenomeAnnotatedDataFrame-method}

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

\section{Coercion to or from other classes}{
  \describe{
    \item{}{

      \code{as(from, "GenomeAnnotatedDataFrame")}:

      Coerce an object of class \code{AnnotatedDataFrame} to a
      \code{GenomeAnnotatedDataFrame}.

    }

    \item{}{
      \code{makeFeatureRanges(object, ...)}:

      Construct a \code{GRanges} instance from a
      \code{GenomeAnnotatedDataFrame}.  Users can specify the genome
      build by passing \code{genome="hg18"} or \code{genome="hg19"}.


    }

  \item{}{
    \code{updateObject(object)}:

    For updating a \code{GenomeAnnotatedDataFrame}

  }

  }
}

\section{Accessors}{
  \describe{
    \item{}{
      \code{chromosome(object)}, \code{chromosome(object) <- value}

      Get or set chromosome.

    }

    \item{}{
      \code{isSnp(object)}:

      Many platforms include polymorphic and nonpolymorphic markers.  isSnp
      evalutes to \code{TRUE} if the marker is polymorphic.

    }


    \item{}{
      \code{position(ojbect)}:

      Physical position in the genome}

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
					       annotationPkg=locusLevelData[["platform"]],
					       genome="hg19")
}
}
\keyword{classes}
