\name{AssayData-methods}
\Rdversion{1.1}
\docType{methods}
\alias{corr}
\alias{flags,AssayData-method}
\alias{nu}
\alias{phi}
\alias{batchNames,AssayData-method}
\alias{batchNames<-,AssayData-method}
\alias{nu,AssayData,character-method}
\alias{phi,AssayData,character-method}

\title{Methods for class AssayData in the oligoClasses package}
\description{

  Batch statistics used for estimating copy number are stored as
AssayData in the 'batchStatistics' slot of the CNSet class. Each element
in the AssayData must have the same number of rows and columns.  Rows
correspond to features and columns correspond to batch.

}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{batchNames}{\code{signature(object = "AssayData")}: ... }
    \item{batchNames<-}{\code{signature(object = "AssayData")}: ... }
    \item{corr}{\code{signature(object = "AssayData", allele = "character")}: ... }
    \item{nu}{\code{signature(object = "AssayData", allele = "character")}: ... }
    \item{phi}{\code{signature(object = "AssayData", allele = "character")}: ... }
%    \item{sigma2}{\code{signature(object = "AssayData", allele = "character")}: ... }
%    \item{tau2}{\code{signature(object = "AssayData", allele = "character")}: ... }
	 }
}

\section{Details}{

\code{lM}:  Extracts entire list of linear model parameters.

\code{corr}: The within-genotype correlation of log2(A) and log2(B) intensities.

\code{nu}: The intercept for the linear model.  The linear model is
fit to the A and B alleles independently.

\code{phi}: The slope for the linear model.  The linear model is fit
independently to the A and B alleles.

%\code{sigma2}: For allele A, sigma2 is calculated as the squared MAD
%of the log2(intensity) for allele 'A' among subjects with genotype AA.
%For allele B, sigma2 is calculated as the squared MAD of the
%log2(intensity) for allele 'B' among subjects with genotype BB.
%sigma2 can be interpreted as a robust estimate of the signal variance.
%
%\code{tau2}: For allele A, tau2 is calculated as the squared MAD of
%the log2(intensity) for allele 'A' among subjects with genotype BB.
%For allele B, tau2 is calculated as the squared MAD of the
%log2(intensity) for allele 'B' among subjects with genotype AA.  tau2
%can be interpeted as a robust estimate of the background variance.

}

\seealso{

	\code{\link{CNSet-class}}

}

\examples{
library(crlmm)
library(Biobase)
data(cnSetExample, package="crlmm")
cnSet <- cnSetExample
isCurrent(cnSet)
assayDataElementNames(batchStatistics(cnSet))
## Accessors for linear model parameters
## -- Included here primarily as a check that accessors are working
## -- Values are all NA until CN estimation is performed using the crlmm package
##
## subsetting
cnSet[1:10, ]
## names of elements in the object
## accessors for parameters
nu(cnSet, "A")[1:10, ]
nu(cnSet, "B")[1:10, ]
phi(cnSet, "A")[1:10, ]
phi(cnSet, "B")[1:10, ]
}

\keyword{classes}
