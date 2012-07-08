\name{SnpSet2-class}
\Rdversion{1.1}
\docType{class}
\alias{SnpSet2-class}
\alias{calls,SnpSet2-method}
\alias{calls<-,SnpSet2,matrix-method}
\alias{confs,SnpSet2-method}
\alias{confs<-,SnpSet2,matrix-method}
\alias{exprs,SnpSet2-method}
\alias{initialize,SnpSet2-method}

\title{Class \code{"SnpSet2"}}

\description{

  A container for genotype calls and confidence scores.  Similar to the
  \code{SnpSet} class in \pkg{Biobase}, but \code{SnpSet2} extends
  \code{gSet} directly whereas \code{SnpSet} extends \code{eSet}.
  Useful properties of \code{gSet} include the \code{genome} slot and
  the \code{GenomeAnnotatedDataFrame}.

}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SnpSet2", assayData, phenoData, featureData, experimentData, annotation, protocolData, call, callProbability, genome, ...)}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
    \item{\code{genome}:}{Object of class \code{"character"} indicating
      the UCSC genome build. Supported builds are 'hg18' and 'hg19'.}
    \item{\code{assayData}:}{Object of class \code{"AssayData"}.}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}.}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}.}
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"}.}
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{gSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "gSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "gSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "gSet", distance 4.
}

\section{Accessors}{

  The argument \code{object} for the following methods is an instance of the
  \code{SnpSet2} class.

  \describe{

    \item{}{
      \code{calls(object)}: \code{calls(object) <- value}:

      Gets or sets the genotype calls. \code{value} can be a
      \code{matrix} or a \code{ff_matrix}.
    }

    \item{}{
      \code{confs(object)}: \code{confs(object) <- value}:

      Gets or sets the genotype confidence scores. \code{value} can be a
      \code{matrix} or a \code{ff_matrix}.
    }

    \item{}{
      \code{snpCall(object)}: \code{snpCallProbability(object) <- value}:

      Gets or sets the genotype confidence scores.
    }
  }
}

\author{
R. Scharpf
}


\seealso{
\code{\linkS4class{SnpSet}}
}
\examples{
showClass("SnpSet2")
new("SnpSet2")
}
\keyword{classes}
