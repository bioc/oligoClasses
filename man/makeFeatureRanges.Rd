\name{makeFeatureGRanges}
\alias{makeFeatureGRanges}

\title{
  Construct a GRanges object from several possible feature-level classes
}
\description{
  Construct a GRanges object from several possible feature-level
  classes. The conversion is useful for subsequent ranged-data queries,
  such as \code{findOverlaps}, \code{countOverlaps}, etc.
}
\usage{
makeFeatureGRanges(object, ...)
}

\arguments{
  \item{object}{
    A \code{gSet}-derived object containing chromosome and physical
  position for the markers on the array.
}

\item{\dots}{
  See the \code{makeFeatureGRanges} method for \code{GenomeAnnotatedDataFrame}.
  }
}

\value{
  A \code{GRanges} object.
}

\author{
R. Scharpf
}

\seealso{
 \code{\link{findOverlaps}}, \code{\linkS4class{GRanges}}, \code{\linkS4class{GenomeAnnotatedDataFrame}}
}
\examples{
	if(require("VanillaICE")){
		library(oligoClasses)
		library(IRanges)
		data(oligoSetExample, package="oligoClasses")
		oligoSet <- oligoSet[chromosome(oligoSet) == 1, ]
		hmmResults <- hmm(oligoSet)

		(frange <- makeFeatureGRanges(oligoSet))
		## which features overlap with the second range in sample NA06993
		olaps <- findOverlaps(hmmResults[["NA06993"]], frange)
		Biobase::featureNames(oligoSet)[subjectHits(olaps)[queryHits(olaps)==2]]

                gr <- makeFeatureGRanges(featureData(oligoSet), genome="hg19")
	}
}
\keyword{manip}

