\name{GenomeAnnotatedDataFrameFrom-methods}
\docType{methods}
\alias{GenomeAnnotatedDataFrameFrom}
\alias{GenomeAnnotatedDataFrameFrom,AssayData-method}
\alias{GenomeAnnotatedDataFrameFrom,ff_or_matrix-method}
\alias{GenomeAnnotatedDataFrameFrom,NULL-method}
\alias{GenomeAnnotatedDataFrameFrom,array-method}
\alias{GenomeAnnotatedDataFrameFrom,list-method}

\title{Methods for Function GenomeAnnotatedDataFrameFrom in Package oligoClasses}
\description{
  \code{GenomeAnnotatedDataFrameFrom} is a convenience for creating
  \code{\linkS4class{GenomeAnnotatedDataFrame}} objects.
}

\section{Methods}{

  Use the method with \code{GenomeAnnotatedDataFrameFrom(object,
    annotationPkg, genome, ...)}; the argument \code{annotationPkg}
    \emph{must} be specified for \code{matrix} and \code{AssayData}
    classes.

\describe{

  \item{\code{signature(object="assayData")}}{

	  This method creates an \code{GenomeAnnotatedDataFrame} using
    feature names and dimensions of an \code{\link{AssayData}} object as
    a template.

  }


  \item{\code{signature(object="matrix")}}{

    This method creates an \code{GenomeAnnotatedDataFrame} using row
    names and dimensions of a \code{\link{matrix}} object as a template.

  }


  \item{\code{signature(object="NULL")}}{

	  This method (called with 'NULL' as the object) creates an
    empty \code{GenomeAnnotatedDataFrame}.

  }

  \item{\code{signature(object="array")}}{

	  This method (called with 'array' as the object) creates a
	  GenomeAnnotatedDataFrame using the first dimension of the
	  array (rows are the number of features).
  }

}}

\author{R Scharpf}

\examples{
	require(Biobase)
	minReqVersion <- "1.0.2"
	require(human370v1cCrlmm)
	if (packageDescription("human370v1cCrlmm", fields='Version') >= minReqVersion){
		x <- matrix(1:25, 5, 5,
			    dimnames=list(c("rs10000092","rs1000055", "rs100016", "rs10003241", "rs10004197"), NULL))
		gd <- GenomeAnnotatedDataFrameFrom(x, annotationPkg="human370v1cCrlmm",
						   genome="hg18")
		pData(gd)
		chromosome(gd)
		position(gd)
	}
}

\keyword{methods}
