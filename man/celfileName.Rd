\name{celfileName}
\alias{celfileName}
\title{
  Extracts complete cel file name from a CNSet object
}
\description{
  Returns the complete cel file (including path) for a CNSet object
}
\usage{
celfileName(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An object of class \code{CNSet}
}
}

\value{
  Character string vector.
}
\author{
  R. Scharpf
}
\note{
  If the CEL files for an experiment are relocated, the \code{datadir}
  should be updated accordingly. See examples.
}

\examples{
  \dontrun{
    if(require(crlmm)){
         data(cnSetExample, package="crlmm")
	 celfileName(cnSetExample)
    }
  }
}
%         datadir=list(dirname="/thumper/ctsa/snpmicroarray/hapmap/raw/affy/1m",
%  n=2L)
%	 fnames <- list.files("/thumper/ctsa/snpmicroarray/hapmap/raw/affy/1m", pattern=".CEL$")
%	 fnames <- fnames[c(grep(sampleNames(cnSetExample)[1], fnames),
%			    grep(sampleNames(cnSetExample)[2], fnames))]
%	 protocolData(cnSetExample)$filename <- fnames
%	 tmp@datadir <- list(dir="/thumper/ctsa/snpmicroarray/hapmap/raw/affy/1m", 2L)
%	 protocolData(tmp) <- protocolData(cnSetExample)
%
%	 data(cnSetExample2)
%	 tmp <- updateObjectFromSlots(cnSetExample2)
%	 tmp@datadir <- list(dir="/thumper/ctsa/snpmicroarray/hapmap/raw/affy/phase3_1m", n=1258L)
%	 protocolData(tmp)$filename <- tmp$celFiles
%	 phenoData(tmp)=phenoData(tmp)[, -4]
%	 cnSetExample2=tmp
%    }
\keyword{manip}
