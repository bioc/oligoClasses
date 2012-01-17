\name{GenomeAnnotatedDataFrameFrom-methods}
\docType{methods}
\alias{GenomeAnnotatedDataFrameFrom}
\alias{GenomeAnnotatedDataFrameFrom,AssayData-method}
\alias{GenomeAnnotatedDataFrameFrom,ff_or_matrix-method}
\alias{GenomeAnnotatedDataFrameFrom,NULL-method}

\title{Methods for Function GenomeAnnotatedDataFrameFrom in Package oligoClasses}
\description{
  \code{GenomeAnnotatedDataFrameFrom} is a convenience for creating
  \code{\linkS4class{GenomeAnnotatedDataFrame}} objects.
}
\section{Methods}{
  Use the method with \code{GenomeAnnotatedDataFrameFrom(object,
    annotationPkg, ...)}; the argument \code{annotationPkg} \emph{must}
  be specified for \code{matrix} and \code{AssayData} classes.
\describe{
  \item{\code{signature(object="assayData")}}{This method creates an
    \code{GenomeAnnotatedDataFrame} using feature names and dimensions of an
    \code{\link{AssayData}} object as a template.}
  \item{\code{signature(object="matrix")}}{This method creates an
    \code{GenomeAnnotatedDataFrame} using row names and dimensions of a
    \code{\link{matrix}} object as a template.}
  \item{\code{signature(object="NULL")}}{This method (called with 'NULL'
    as the object) creates an empty \code{GenomeAnnotatedDataFrame}.}
}}

\author{R Scharpf}

\keyword{methods}