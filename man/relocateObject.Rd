\name{relocateObject}
\alias{relocateObject}
\title{ Safely moves all files saved in PATH to a new directory by updating the pointers inside of the ff objects.}
\description{
	See details
}
\usage{
relocateObject(object, ...)
}

\arguments{
  \item{object}{
	  An object of class \code{CNSet}.
  }
  \item{\dots}{
	  Not used.
  }
}
\details{

When crlmm is used in conjunction with the ff package, large files are
saved as ff_matrix of ffdf objects to avoid excessive RAM.  The R
object in RAM contains pointers to the location on disk.  If the files
on disk are moved, the data can no longer be accessed through the R
object in memory.  However, if the pointers contained in the R object
are modified to point to the new location, the data can be accessed as
usual.  The relocate function edits the pointers of all ff objects
embedded within the CNSet class. In particular, all elements stored in
slot assayData and all elements stored in slot batchStatistics.

}
\value{
	An object of class \code{CNSet}.
}

\author{
R. Scharpf and L. Mireless
}
\seealso{
	\code{\link{CNSet-class}}
}
\examples{
	\dontrun{
		library(ff)
		library(crlmm)
		cdfName <- "genomewidesnp6"
		pathToCels <- "/thumper/ctsa/snpmicroarray/hapmap/raw/affy/1m"
		outdir <- file.path(getwd(), "ff")
		if(!file.exists(outdir)) stop("Please specify a valid directory for storing output")
		if(!file.exists(pathToCels)) stop("Please specify the correct path to the CEL files")
		celFiles <- list.celfiles(pathToCels, full.names=TRUE, pattern=".CEL")[1:6]
		batch <- as.factor(substr(basename(celFiles), 13, 13))
		ldPath(outdir)
		ocProbesets(100e3)
		ocSamples(200)
		gtSet(celFiles, cdfName, batch=batch)
		A(gtSet)[1:5, ]
		pattern(A(gtSet))
		filename(A(gtSet))

		outdir2 <- file.path(getwd(), "testing")
		relocateObject(gtSet, from=outdir,
			       to=outdir2)
		## note the pointers have changed
		pattern(A(gtSet))
		filename(A(gtSet))
	}
}
\keyword{manip}

