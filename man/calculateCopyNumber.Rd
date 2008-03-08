\name{calculateCopyNumber}
\alias{calculateCopyNumber}
\title{Calculates copy number as the average of SNP-level summaries from
  the A and B alleles}

\description{ Calculates copy number as the average of SNP-level
summaries from the A and B alleles from an object of class \code{SnpCallSetPlus} }

\usage{
calculateCopyNumber(object, ...)
}

\arguments{
  \item{object}{an instance of \code{SnpCallSetPlus}}
  \item{\dots}{Not implemented}
}

\details{

  The average A (and B) allele intensities for each SNP is centered by
  the median intensity across samples (should have a reasonable (e.g., >
  10) samples).  Addtionally, we assume that the median intensity across
  all samples for autosomes should be the same for homozygous genotypes
  and heterozygous geneotypes.  The resulting intensities are then
  recentered at 2 copies (for autosomes), 1 for the male X chromosome, 2
  for the female X, etc.  We assume that for a given SNP, the median
  copy number is normal. If this is not likely to be the case,
  alternative approaches should be explored.

  }

\value{

  A matrix of copy number intensities.  

}

\author{R. Scharpf}

\seealso{\code{\link{SnpCallSetPlus-class}}}
\keyword{methods}

