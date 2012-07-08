\name{BeadStudioSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{BeadStudioSetList-class}
\alias{oligoSetList-class}
\title{BeadStudioSetList class}

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

\author{
R. Scharpf
}
\seealso{
  See supporting packages for methods defined for the class.
}

\keyword{classes}