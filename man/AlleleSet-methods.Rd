\name{getA}
\alias{getA}
\alias{getM}
\alias{getA,SnpQSet-method}
\alias{getM,SnpQSet-method}
\alias{getA,SnpCnvQSet-method}
\alias{getM,SnpCnvQSet-method}
\alias{getM,TilingFeatureSet2-method}
\alias{getA,TilingFeatureSet2-method}

\alias{A}
\alias{A,AlleleSet-method}
\alias{A<-}
\alias{A<-,AlleleSet-method}
\alias{A<-,AlleleSet,matrix-method}
\alias{B}
\alias{B,AlleleSet-method}
\alias{B<-}
\alias{B<-,AlleleSet-method}
\alias{B<-,AlleleSet,matrix-method}
\alias{open,AlleleSet-method}
\alias{close,AlleleSet-method}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute average log-intensities / log-ratios}
\description{
  Methods to compute average log-intensities and log-ratios across
  alleles, within strand.
}
\usage{
getA(object)
getM(object)
A(object, ...)
B(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{SnpQSet}, \code{SnpCnvQSet} or
    \code{TilingFeatureSet2} object.}
%  \item{con}{\code{AlleleSet} or \code{AlleleSet} extension.}
  \item{...}{arguments to be passed to \code{allele} - 'sense' and
    'antisense' are valid values if the array is pre-SNP_5.0}
}

\details{

  For SNP data, SNPRMA summarizes the SNP information into 4 quantities (log2-scale):
  \itemize{
    \item{antisenseThetaA}{antisense allele A. (Not applicable for Affymetrix 5.0 and 6.0 platforms.)}
    \item{antisenseThetaB}{antisense allele B. (Not applicable for Affymetrix 5.0 and 6.0 platforms.)}
    \item{senseThetaA}{sense allele A.  (Not applicable for Affymetrix 5.0 and 6.0 platforms.)}
    \item{senseThataB}{sense allele B.  (Not applicable for Affymetrix 5.0 and 6.0 platforms.)}
    \item{alleleA}{Affymetrix 5.0 and 6.0 platforms}
    \item{alleleB}{Affymetrix 5.0 and 6.0 platforms}
  }

  The average log-intensities are given by:
  \code{(antisenseThetaA+antisenseThetaB)/2} and
  \code{(senseThetaA+senseThetaB)/2}.

  The average log-ratios are given by:
  \code{antisenseThetaA-antisenseThetaB} and
  \code{senseThetaA-senseThetaB}.

  For Tiling data, \code{getM} and \code{getA} return the log-ratio and
  average log-intensities computed across channels:
  \code{M = log2(channel1)-log2(channel2)}
  \code{A = (log2(channel1)+log2(channel2))/2}

  When large data support is enabled with the ff package, the
  AssayData elements of an \code{AlleleSet} object can be
  \code{ff_matrix} or \code{ffdf}, in which case pointers to the ff
  object are stored in the assay data.  The functions \code{open} and
  \code{close} can be used to open or close the connection,
  respectively.

}
\value{
  A 3-dimensional array (SNP's x Samples x Strand) with the requested
  measure, when the input SNP data (50K, 250K).

  A 2-dimensional array (SNP's x Samples), when the input is from SNP
  5.0 and SNP 6.0 arrays.

  A 2-dimensional array if the input is from Tiling arrays.
}

\seealso{\code{\link[oligo]{snprma}}}
\keyword{manip}

