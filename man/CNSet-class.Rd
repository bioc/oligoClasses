\name{CNSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CNSet-class}
\alias{CNSet}
\alias{[,CNSet-method}
\alias{[,CNSet,ANY-method}
\alias{A,CNSet-method}
\alias{A<-,CNSet-method}
\alias{allele,CNSet-method}
\alias{B,CNSet-method}
\alias{B<-,CNSet-method}
\alias{batch,CNSet-method}
\alias{batchNames,CNSet-method}
\alias{batchNames<-,CNSet-method}
\alias{batchStatistics,CNSet-method}
\alias{batchStatistics<-,CNSet,AssayData-method}
\alias{calls,CNSet-method}
\alias{calls<-,CNSet,matrix-method}
\alias{close,CNSet-method}
\alias{coerce,CNSetLM,CNSet-method}
\alias{coerce,CNSet,CopyNumberSet-method}
\alias{coerce,CNSet,oligoSnpSet}
\alias{coerce,CNSet,oligoSnpSet-method}
\alias{confs,CNSet-method}
\alias{confs<-,CNSet,matrix-method}
\alias{corr,CNSet,character-method}
\alias{flags,CNSet-method}
\alias{initialize,CNSet-method}
\alias{initialize,CNSetLM-method}
%\alias{lM,CNSet-method}
\alias{nu,CNSet,character-method}
\alias{open,CNSet-method}
\alias{phi,CNSet,character-method}
\alias{show,CNSet-method}
\alias{snpCallProbability,CNSet-method}
\alias{sigma2,CNSet,character-method}
\alias{tau2,CNSet,character-method}
\alias{updateObject,CNSet-method}

\title{Class "CNSet"}
\description{

	CNSet is a container for intermediate data and parameters
	 pertaining to allele-specific copy number estimation. Methods
	 for CNSet objects, including accessors for linear model
	 parameters and allele-specific copy number are included here.

}
\section{Objects from the Class}{

		 An object from the class is not generally intended to
		 be initialized by the user, but returned by the
		 \code{genotype} function in the \code{crlmm} package.

		 The following creates a very basic \code{CNSet} with
		 \code{assayData} containing the required elements.

	 \code{new(CNSet, alleleA=new("matrix"), alleleB=new("matrix"), call=new("matrix"),
		   callProbability=new("matrix"), batch=new("factor"))}

}
\section{Slots}{
  \describe{
    \item{\code{batch}:}{Object of class \code{"factor"} ~~ }
    \item{\code{batchStatistics}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class
      \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{datadir}:}{Object of class \code{"list"}~~}
    \item{\code{mixtureParams}:}{Object of class \code{"matrix"}~~}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}

\section{Methods}{

  The argument \code{object} for the following methods is a \code{CNSet}.

  \describe{
    \item{}{\code{object[i, j]}: subset the \code{CNSet} object by
    markers (i) and/or samples (j). }

    \item{}{\code{A(objet)}: accessor for the normalized intensities of
    allele A}

  \item{}{\code{A(object) <- value}: replace intensities for the A
    allele intensities by \code{value}. The object \code{value} must be
    a \code{matrix}, \code{ff_matrix}, or \code{ffdf}.}


    \item{}{\code{allele(object, allele)}: accessor for the normalized
    intensities for the A or B allele. The argument for \code{allele}
    must be either 'A' or 'B'}

  \item{}{\code{B(objet)}: accessor for the normalized intensities of
    allele B}

  \item{}{\code{B(object) <- value}: replace intensities for the B
    allele intensities by \code{value}. The object \code{value} must be
    a \code{matrix}, \code{ff_matrix}, or \code{ffdf}.}

  \item{}{\code{batch(object)}: vector of batch labels for each sample.}

  \item{}{\code{batchNames(object)}: the unique batch names}

  %\item{batchNames<-}{\code{signature(object = "CNSet")}: ... }

  \item{}{\code{batchNames(object) <- value}: relabel the batches}

  \item{}{\code{calls(object)}: accessor for genotype calls coded as 1
    (AA), 2 (AB), or 3 (BB).  Nonpolymorphic markers are \code{NA}.}

    \item{}{\code{confs(object)}: accessor for the genotype confidence scores.}

  \item{}{\code{close(object)}: close any open file connections to
  \code{ff} objects stored in the \code{CNSet} object.}

   %\item{coerce}{\code{signature(from="CNSet")}: ... }
   \item{}{\code{as(object, "oligoSnpSet")}: coerce a \code{CNSet}
   object to an object of class \code{oligoSnpSet} -- a container for
   the total copy number and genotype calls.}

%    \item{corr}{\code{signature(object = "CNSet", allele =
%	"character")}: ... }

    \item{}{\code{corr(object)}: the correlation of the A and B
    intensities within each genotype.}

  %\item{flags}{\code{signature(object="CNSet")}: SNP flags }
  \item{}{\code{flags(object)}:  flags to indicate possible problems with
  the copy number estimation.  Not fully implemented at this point.}

   \item{}{\code{new("CNSet")}: instantiating a \code{CNSet} object.}

   \item{}{\code{nu(object, allele)}: accessor for the intercept
   (background) for the A and B alleles.  The value of \code{allele}
   must be 'A' or 'B'.}

  \item{}{\code{open(object)} open file connections for all \code{ff}
    objects stored in the \code{CNSet} object.}

   \item{}{\code{nu(object, allele)}: accessor for the slope for the A
    and B alleles.  The value of \code{allele} must be 'A' or 'B'.}

   \item{}{\code{sigma2(object, allele)}: accessor for the within
   genotype variance}

   \item{}{\code{tau2(object, allele)}: accessor for background variance}
 }
 }

\author{
R. Scharpf
}

\examples{
	new("CNSet")
}
\keyword{classes}
