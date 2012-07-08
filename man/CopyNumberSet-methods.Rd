\name{CopyNumberSet-methods}
\docType{methods}
\alias{CopyNumberSet-methods}
\alias{copyNumber}
\alias{copyNumber<-}
\alias{cnConfidence}
\alias{cnConfidence<-}
\title{Methods for class CopyNumberSet.}
\description{
	Accessors and CopyNumberSet
}

\usage{

copyNumber(object, ...)
cnConfidence(object)
copyNumber(object) <- value
cnConfidence(object) <- value

}

\arguments{

\item{object}{\code{CopyNumberSet} object or derived class}
\item{\dots}{Ignored for \code{CopyNumberSet} and \code{oligoSnpSet}. }
\item{value}{matrix}
}

\value{

	\code{copyNumber} returns a matrix of copy number estimates or
	relative copy number estimates.  Since the copy number estimates
	are stored as integers (copy number * 100), the matrix returned
	by the \code{copyNumber} accessor will need to be divided by a
	factor of 100 to transform the measurements back to the original
	copy number scale.

	\code{cnConfidence} returns a matrix of confidence scores for
	the copy number estimates.  These are also represented as
	integers and will require a back-transformation to the original
	scale.

}
\keyword{methods}
\keyword{manip}

\examples{
	library(Biobase)
	data(locusLevelData)
	path <- system.file("extdata", package="oligoClasses")
	fd <- readRDS(file.path(path, "genomeAnnotatedDataFrameExample.rds"))
	## the following command creates an 'oligoSnpSet' object, storing
	## an integer representation of the log2 copy number in the 'copyNumber' element
	## of the assayData.  Genotype calls and genotype confidence scores are also stored
	## in the assayData.
	oligoSet <- new("oligoSnpSet",
			copyNumber=integerMatrix(log2(locusLevelData[["copynumber"]]/100), 100),
			call=locusLevelData[["genotypes"]],
			callProbability=integerMatrix(locusLevelData[["crlmmConfidence"]], 1),
			annotation=locusLevelData[["platform"]],
			featureData=fd,
			genome="hg19")

	## There are several accessors for the oligoSnpSet class.
	icn <- copyNumber(oligoSet)
	range(icn) ## integer scale
	lcn <- icn/100
	range(lcn) ## log2 copy number

	## confidence scores for the genotypes are also represented on an integer scale
	ipr <- snpCallProbability(oligoSet)
	range(ipr) ## integer scale

	## for genotype confidence scores, the helper function i2p
	## converts back to a probability scale
	pr <- i2p(ipr)
	range(pr)

	## The helper function confs is a shortcut, extracting the
	## integer-based confidence scores and transforming to the
	## probability scale
	pr2 <- confs(oligoSet)
	all.equal(pr, pr2)

	## To extract information on the annotation of the SNPs, one can use
	position(oligoSet)
	chromosome(oligoSet)
	## the position and chromosome coordinates were extracted from build hg19
	genomeBuild(oligoSet)
}

