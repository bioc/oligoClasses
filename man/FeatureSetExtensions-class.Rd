\name{FeatureSet-class}
\docType{class}
\alias{show,FeatureSet-method}
\alias{FeatureSet-class}
\alias{initialize,FeatureSet-method}
\alias{GeneFeatureSet-class}
\alias{ExpressionFeatureSet-class}
\alias{ExonFeatureSet-class}
\alias{SnpFeatureSet-class}
\alias{SnpCnvFeatureSet-class}
\alias{TilingFeatureSet-class}

\title{"FeatureSet" and "FeatureSet" Extensions}
\description{Classes to store data from Expression/Exon/SNP/Tiling arrays at the feature level.}
\section{Objects from the Class}{
  The FeatureSet class is VIRTUAL. Therefore users are not able to
  create instances of such class.
  
  Objects for FeatureSet-like classes can be created by calls of the form:
  \code{new(CLASSNAME, assayData, manufacturer, platform, exprs,
    phenoData, featureData, experimentData, annotation, ...)}.
  But the preferred way is using parsers like
  \code{\link[oligo]{read.celfiles}} and \code{\link[oligo]{read.xysfiles}}.
}
\section{Slots}{
  \describe{
    \item{\code{manufacturer}:}{Object of class \code{"character"}  }
    \item{\code{assayData}:}{Object of class \code{"AssayData"}  }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}  }
  }
}

\section{Methods}{
 \describe{
    \item{show}{\code{signature(.Object = "FeatureSet")}: show object contents }
    }
}
\author{Benilton Carvalho}

\seealso{\code{\link[Biobase:class.eSet]{eSet}}, \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, \code{\link[Biobase:class.Versioned]{Versioned}}}

\examples{
set.seed(1)
tmp <- 2^matrix(rnorm(100), ncol=4)
rownames(tmp) <- 1:25
colnames(tmp) <- paste("sample", 1:4, sep="")
efs <- new("ExpressionFeatureSet", exprs=tmp)
}
\keyword{classes}
