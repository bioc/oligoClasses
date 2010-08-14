\name{LinearModelParameter-class}
\Rdversion{1.1}
\docType{class}
\alias{corr}
\alias{nu}
\alias{phi}
\alias{sigma2}
\alias{tau2}
\alias{lM}
\alias{lM<-}
\alias{LinearModelParameter-class}
\alias{batchNames,LinearModelParameter-method}
\alias{batchNames<-,LinearModelParameter-method}
\alias{corr,LinearModelParameter,character-method}
\alias{lM<-,CNSet,LinearModelParameter-method}
\alias{nu,LinearModelParameter,character-method}
\alias{phi,LinearModelParameter,character-method}
\alias{sigma2,LinearModelParameter,character-method}
\alias{tau2,LinearModelParameter,character-method}

\title{Class "LinearModelParameter"}
\description{

	A container for storing the estimated parameters for copy
	number. This container is a union of 'AssayData', 'list', and
	'environment' classes.  Each element must have the same number
	of rows and columns.  Rows correspond to features and columns
	correspond to batch.

}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{batchNames}{\code{signature(object = "LinearModelParameter")}: ... }
    \item{batchNames<-}{\code{signature(object = "LinearModelParameter")}: ... }
    \item{corr}{\code{signature(object = "LinearModelParameter", allele = "character")}: ... }
    \item{lM<-}{\code{signature(object = "CNSet", value = "LinearModelParameter")}: ... }
    \item{nu}{\code{signature(object = "LinearModelParameter", allele = "character")}: ... }
    \item{phi}{\code{signature(object = "LinearModelParameter", allele = "character")}: ... }
    \item{sigma2}{\code{signature(object = "LinearModelParameter", allele = "character")}: ... }
    \item{tau2}{\code{signature(object = "LinearModelParameter", allele = "character")}: ... }
	 }
}

\seealso{

	\code{\link{CNSet-class}}	

}

\examples{

x <- matrix(runif(250*96*2, 0, 2), 250, 96*2)
test1 <- new("CNSet", alleleA=x, alleleB=x, call=x, callProbability=x,
	     batch=as.factor(rep(letters[1:2], each=96)))
## Accessors for linear model parameters
## -- Included here primarily as a check that accessors are working
## -- Values are all NA until CN estimation is performed using the crlmm package
##
## subsetting
test1[1:10, 1:5]
## names of elements in the object
ls(lM(test1))
## accessors for parameters
nu(test1, "A")[1:10, ]
nu(test1, "B")[1:10, ]
phi(test1, "A")[1:10, ]
phi(test1, "B")[1:10, ]
sigma2(test1, "A")[1:10, ]
sigma2(test1, "B")[1:10, ]
tau2(test1, "A")[1:10, ]
tau2(test1, "B")[1:10, ]
corr(test1, "AA")[1:10, ]
corr(test1, "AB")[1:10, ]
corr(test1, "BB")[1:10, ]


}

\keyword{classes}
