\name{gSet-class}
\Rdversion{1.1}
\docType{class}
\alias{gSet-class}
\alias{gSet}
\alias{checkOrder,gSet-method}
\alias{db,gSet-method}
\alias{genomeBuild,gSet-method}
\alias{genomeBuild<-,gSet,character-method}
\alias{getArm}
\alias{getArm,gSet-method}
\alias{initialize,gSet-method}
\alias{makeFeatureGRanges,gSet-method}
\alias{show,gSet-method}

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
    \item{\code{protocolData}:}{Object of class
      \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{genome}:}{Object of class \code{"character"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  The \code{object} for the below methods is a class that extends the
  virtual class \code{gSet}.
  \describe{
    \item{}{\code{checkOrder(object)}: checks that the object is ordered
      by chromosome and physical position. Returns \code{logical}.}
    \item{}{\code{chromosome(object)}: accessor for chromosome in the
      \code{GenomeAnnotatedDataFrame} slot.}
    \item{}{\code{chromosome(object) <- value}: replacement method for chromosome in the
      \code{GenomeAnnotatedDataFrame} slot. \code{value} must be an
      \code{integer} vector.}
    \item{}{
      \code{db(object)}: database connection
    }
    \item{}{

       \code{genomeBuild(object)}, \code{genomeBuild(object) <- value}:

        Get or set the UCSC genome build. Supported builds are hg18 and hg19.

      }

      \item{}{

	\code{getArm(object)}:  Character vector indicating the chromosomal arm for
	each marker in \code{object}.

      }


    \item{}{\code{isSnp(object)}: whether the marker is
      polymorphic. Returns a \code{logical} vector.}

    \item{}{

      \code{makeFeatureGRanges(object)}:  Construct an instance of the
      \code{GRanges} class from a \code{GenomeAnnotatedDataFrame}.

    }

    \item{}{\code{position(object)}: \code{integer} vector of the
      genomic position}
    \item{}{
      \code{show(object)}:

      Print a concise summary of \code{object}.

    }
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
