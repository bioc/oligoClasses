\name{DBPDInfo-class}
\docType{class}
\alias{initialize,DBPDInfo-method}
\alias{annotation,DBPDInfo-method}
\alias{DBPDInfo-class}
\alias{AffySNPPDInfo-class}
\alias{AffySNPCNVPDInfo-class}
\alias{show,DBPDInfo-method}
\alias{AffyExpressionPDInfo-class}
\alias{AffyTilingPDInfo-class}
\alias{AffyGenePDInfo-class}
\alias{AffyExonPDInfo-class}
\alias{AffySTPDInfo-class}
\alias{SNPPDInfo-class}
\alias{SNPCNVPDInfo-class}
\alias{NgsTilingPDInfo-class}
\alias{TilingPDInfo-class}
\alias{ExpressionPDInfo-class}
\alias{NgsExpressionPDInfo-class}

\alias{DBPDInfo}
\alias{AffySNPPDInfo}
\alias{AffySNPCNVPDInfo}
\alias{show,DBPDInfo-method}
\alias{AffyExpressionPDInfo}
\alias{AffyTilingPDInfo}
\alias{AffyGenePDInfo}
\alias{AffyExonPDInfo}
\alias{AffySTPDInfo}
\alias{SNPPDInfo}
\alias{SNPCNVPDInfo}
\alias{NgsTilingPDInfo}
\alias{TilingPDInfo}
\alias{ExpressionPDInfo}
\alias{NgsExpressionPDInfo}

\title{Class "DBPDInfo"}
\description{A class for Platform Design Information objects, stored
  using a database approach}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DBPDInfo", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{getdb}:}{Object of class \code{"function"}  }
    \item{\code{tableInfo}:}{Object of class \code{"data.frame"}  }
    \item{\code{manufacturer}:}{Object of class \code{"character"}  }
    \item{\code{genomebuild}:}{Object of class \code{"character"}  }
    \item{\code{geometry}:}{Object of class \code{"integer"} with length
    2 (rows x columns)}
  }
}
\section{Methods}{
  \describe{
    \item{annotation}{string describing annotation package associated to
      object}
    }
}
\keyword{classes}
