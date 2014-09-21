\name{BeadStudioSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{BeadStudioSetList-class}
\alias{BafLrrSetList-class}
\alias{oligoSetList-class}
\alias{[[,BeadStudioSetList,ANY,ANY-method}
\alias{[[,BafLrrSetList,ANY,ANY-method}
\alias{[[<-,BafLrrSetList,ANY,ANY,BafLrrSet-method}
\alias{baf,BeadStudioSetList-method}
\alias{baf,oligoSetList-method}
\alias{calls,oligoSetList-method}
\alias{clone}
\alias{clone,BafLrrSetList-method}
\alias{copyNumber,oligoSetList-method}
\alias{initialize,BeadStudioSetList-method}
\alias{lrr,BeadStudioSetList-method}
\alias{lrr<-,BafLrrSetList,matrix-method}
\alias{updateObject,BeadStudioSetList-method}

\title{List classes with assay data listed by chromosome}

\description{
  Container for log R ratios and B allele frequencies stored by chromosome.
}


\section{Slots}{
  \describe{
    \item{\code{assayDataList}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureDataList}:}{Object of class \code{"list"} ~~ }
    \item{\code{chromosome}:}{Object of class \code{"integer"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{genome}:}{Object of class \code{"character"} indicating
  the genome build. Valid entries are "hg18" and "hg19". }
  }
}

\section{Methods defined for the class}{
  \describe{
    \item{}{\code{clone(object, id, prefix="",...)}

      Performs a deep copy of the ff objects in the assay data elements
  of \code{object}.  A new object of the same class will be
  instantiated.  The ff objects in the instantiated object will point to
  ff files on disk with prefix given by the argument \code{prefix}.

  A use-case for such a function is that one may want to perform wave
  correction on the log R ratios in \code{object}, but keep a copy of
  the original unadjusted log R ratios.  If \code{object} is not copied
  using \code{clone} prior to wave correction, the log R ratios will be
  updated on disk and the original, unadjusted log R ratios will no
  longer be available.

  }
}
}

\author{
R. Scharpf
}
\seealso{
  See supporting packages for methods defined for the class.
}

\keyword{classes}