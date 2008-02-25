\name{calculateCopyNumber}
\alias{calculateCopyNumber}
\title{Calculates copy number as the average of SNP-level summaries from
  the A and B alleles}

\description{ Calculates copy number as the average of SNP-level
summaries from the A and B alleles from an object of class \code{SnpCallSetPlus} }

\usage{
calculateCopyNumber(object)
}

\arguments{
  \item{object}{an instance of \code{SnpCallSetPlus}}
}

\value{

  A matrix of copy number intensities.  The average of the A and B
  allele intensities for a sample is centered at log2(1) -- this assumes
  that the average copy number in a sample is two.  

}

\author{R. Scharpf}

\seealso{\code{\link{SnpCallSetPlus-class}}}
\keyword{methods}

