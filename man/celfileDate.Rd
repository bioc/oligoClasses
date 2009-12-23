\name{celfileDate}
\alias{celfileDate}
\title{ Cel file dates}
\description{ Parses cel file dates from the header of .CEL files for the Affymetrix platform}
\usage{
	celfileDate(filename)
}
\arguments{
  \item{filename}{
  Name of cel file
}
}
\value{
	character string
}
\author{H. Jaffee}
\examples{
require(hapmapsnp6)
path <- system.file("celFiles", package="hapmapsnp6")
celfiles <- list.celfiles(path, full.names=TRUE)
dts <- sapply(celfiles, celfileDate)
}
\keyword{manip}

