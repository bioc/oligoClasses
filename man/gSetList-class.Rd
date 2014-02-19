\name{gSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{gSetList-class}
\alias{$,gSetList-method}
\alias{$<-,gSetList-method}
\alias{[,gSetList-method}
\alias{[,gSetList,ANY-method}
\alias{[[<-,gSetList,ANY,ANY,BafLrrSet-method}
\alias{annotation,gSetList-method}
\alias{assayDataList,gSetList-method}
\alias{chromosome,gSetList-method}
\alias{coerce,gSetList,list-method}
\alias{dims,gSetList-method}
\alias{elementLengths,gSetList-method}
\alias{featureDataList,gSetList-method}
\alias{genomeBuild,gSetList-method}
\alias{genomeBuild<-,gSetList,character-method}
\alias{initialize,gSetList-method}
\alias{length,gSetList-method}
\alias{makeFeatureGRanges,gSetList-method}
\alias{position,gSetList-method}
\alias{sampleNames,gSetList-method}
\alias{sampleNames<-,gSetList,character-method}
\alias{show,gSetList-method}


\title{Virtual Class for Lists of eSets}
\description{
  Virtual Class for Lists of eSets.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{assayDataList}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{featureDataList}:}{Object of class \code{"list"} ~~ }
    \item{\code{chromosome}:}{Object of class \code{"vector"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{genome}:}{Object of class \code{"character"} ~~ }
  }
}

\section{Accessors}{

  \code{object}  is an instance of a \code{gSetList}-derived class.

  \describe{
    \item{}{
      \code{annotation(object)}:

      character string indicating the package used to provide annotation
      for the features on the array.
    }

    \item{}{
      \code{chromosome(object)}:

      Returns the chromosome corresponding to each element in the
      \code{gSetList} object

    }

    \item{}{
      \code{elementLengths(object)}: Returns the number of rows for each
      list of assays.  In most \code{gSetList}-derived classes, the
      assays are organized by chromosome and \code{elementLengths}
      returns the number of markers for each chromosome.
    }


    \item{}{
      \code{genomeBuild(object)}, \code{genomeBuild(object) <- value}:

      Get or set the UCSC genome build. Supported builds are hg18 and hg19.

    }

  }

}

\section{Coercion}{
  \code{object}  is an instance of a \code{gSetList}-derived class.

  \describe{
    \item{}{
      \code{makeFeatureGRanges(object, ...)}:

      Create a GRanges object for the featureData.  The featureData is
      stored as a list.  This method stacks the featureData from each
      list element.  Metadata columns in the GRanges object include
      physical position ('position'), a SNP indicator ('isSnp'), and the
      chromosome. The genome build is extracted from \code{object} using
      the method \code{genomeBuild}.

      }

  }

}

\author{
R. Scharpf
}

\seealso{
  \code{\linkS4class{oligoSetList}}, \code{\linkS4class{BeadStudioSetList}}
}
\examples{
showClass("gSetList")
}
\keyword{classes}
