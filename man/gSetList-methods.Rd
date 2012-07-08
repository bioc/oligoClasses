\name{gSetList-methods}
\docType{methods}
\alias{annotation,gSetList-method}
\alias{genomeBuild,gSetList-method}
\alias{genomeBuild<-,gSetList,character-method}
\alias{initialize,gSetList-method}

\title{Methods for gSetList objects}

\description{Methods for gSetList objects}

\section{Accessors}{

  \code{object}  is an instance of a \code{gSetList}-derived class.

  \describe{
    \item{}{
      \code{annotation(object)}:

      character string indicating the package used to provide annotation
      for the features on the array.
    }

    \item{}{
      \code{genomeBuild(object)}, \code{genomeBuild(object) <- value}:

      Get or set the UCSC genome build. Supported builds are hg18 and hg19.

    }

  }

}

\keyword{methods}
