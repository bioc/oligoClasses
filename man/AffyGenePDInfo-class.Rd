\name{AffyGenePDInfo-class}
\docType{class}
\alias{AffyGenePDInfo-class}

\title{Class AffyGenePDInfo}
\description{Short description}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AffyGenePDInfo", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{getdb}:}{Object of class \code{"function"} }
    \item{\code{tableInfo}:}{Object of class \code{"data.frame"} }
    \item{\code{geometry}:}{Object of class \code{"integer"} }
    \item{\code{manufacturer}:}{Object of class \code{"character"} }
    \item{\code{genomebuild}:}{Object of class \code{"character"} }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{AffyExpressionPDInfo}"}, directly.
  Class \code{"\linkS4class{ExpressionPDInfo}"}, by class "AffyExpressionPDInfo", distance 2.
  Class \code{"\linkS4class{DBPDInfo}"}, by class "AffyExpressionPDInfo", distance 3.
  Class \code{"\linkS4class{PDInfo}"}, by class "AffyExpressionPDInfo", distance 4.
}
\section{Methods}{
No methods defined with class "AffyGenePDInfo" in the signature.
}
\author{RS}
\seealso{
	\code{\link[oligoClasses:class.DBPDInfo]{DBPDInfo}}
}
	
\examples{
showClass("AffyGenePDInfo")
}
\keyword{classes}
