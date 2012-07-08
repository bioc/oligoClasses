\name{makeFeatureRanges}
\alias{makeFeatureRanges}

\title{
  Construct a GRanges object from several possible feature-level classes
}
\description{
  Construct a GRanges object from several possible feature-level
  classes. The conversion is useful for subsequent ranged-data queries,
  such as \code{findOverlaps}, \code{countOverlaps}, etc.
}
\usage{
makeFeatureRanges(object)
}

\arguments{
  \item{object}{
    A \code{gSet}-derived object containing chromosome and physical
  position for the markers on the array.
  }
}

\value{
  A \code{GRanges} object.
}

\author{
R. Scharpf
}

\seealso{
 \code{\link{findOverlaps}}, \code{\linkS4class{GRanges}}
}
\examples{
	if(require("VanillaICE")){
		library(oligoClasses)
		library(IRanges)
		data(oligoSetExample, package="oligoClasses")
		oligoSet <- oligoSet[chromosome(oligoSet) == 1, ]
		hmmResults <- hmm(oligoSet)

		(frange <- makeFeatureRanges(oligoSet))
		## which features overlap with the second range in sample NA06993
		olaps <- findOverlaps(hmmResults[["NA06993"]], frange)
		Biobase::featureNames(oligoSet)[subjectHits(olaps)[queryHits(olaps)==2]]
	}
}
\keyword{manip}

