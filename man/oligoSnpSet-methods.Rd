\name{oligoSnpSet-methods}
\docType{methods}
\alias{oligoSnpSet-class}
\alias{initialize,oligoSnpSet-method}
\alias{baf,oligoSnpSet-method}
\alias{baf<-,oligoSnpSet-method}
\alias{coerce,oligoSnpSet,data.frame-method}
\alias{copyNumber,oligoSnpSet-method}
\alias{copyNumber<-,oligoSnpSet,matrix-method}
\alias{calls,oligoSnpSet-method}
\alias{calls<-,oligoSnpSet,matrix-method}
\alias{cnConfidence,oligoSnpSet-method}
\alias{cnConfidence<-,oligoSnpSet,matrix-method}
\alias{callsConfidence,oligoSnpSet-method}
\alias{callsConfidence<-,oligoSnpSet,matrix-method}
\alias{position<-,oligoSnpSet,integer-method}
\alias{updateObject,oligoSnpSet-method}
\title{Methods for oligoSnpSet class}
\description{Methods for oligoSnpSet class}

\section{Methods}{

  In the following code, \code{object} is an instance of the
  \code{oligoSnpSet} class.

  \describe{

    \item{}{
      \code{new("oligoSnpSet", ...)}:
      Instantiates an object of class \code{oligoSnpSet}.  The assayData
      elements of the \code{oligoSnpSet} class can include matrices of
      genotype calls, confidence scores for the genotype calls, B allele
      frequencies, absolute or relative copy number, and confidence
      scores for the copy number estimates. Each matrix should be
      coerced to an integer scale prior to assignment to the
      \code{oligoSnpSet} object.  Validity methods defined for the class
      will fail if the matrices are not integers. See examples for
      additional details.
    }

    \item{}{
      \code{baf(object)}:
      Accessor for integer representation of the B allele frequencies.
      The value returned by this method can be divided by 1000 to obtain
      B allele frequencies on the original [0,1] scale.
    }

    \item{}{
      \code{baf(object) <- value}:
      Assign an integer representation of the B allele frequencies to
      the 'baf' element of the assayData slot.  \code{value} must be a
      matrix of integers. See the examples for help converting BAFs to a
      matrix of integers.
    }
  }

}

%\examples{}


\keyword{methods}
