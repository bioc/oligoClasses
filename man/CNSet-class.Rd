\name{CNSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CNSet-class}
\alias{[,CNSet-method}
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
\alias{close,CNSet-method}
\alias{coerce,CNSetLM,CNSet-method}
\alias{corr,CNSet,character-method}
\alias{flags,CNSet-method}
\alias{initialize,CNSet-method}
\alias{initialize,CNSetLM-method}
%\alias{lM,CNSet-method}
\alias{nu,CNSet,character-method}
\alias{open,CNSet-method}
\alias{phi,CNSet,character-method}
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
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SnpSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "SnpSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "SnpSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "SnpSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "CNSet")}: ... }
    \item{A}{\code{signature(object = "CNSet")}: ... }
    \item{A<-}{\code{signature(object = "CNSet")}: ... }
    \item{allele}{\code{signature(object = "CNSet")}: ... }
    \item{B}{\code{signature(object = "CNSet")}: ... }
    \item{B<-}{\code{signature(object = "CNSet")}: ... }
    \item{batch}{\code{signature(object = "CNSet")}: ... }
    \item{batchNames}{\code{signature(object = "CNSet")}: ... }
    \item{batchNames<-}{\code{signature(object = "CNSet")}: ... }
    \item{close}{\code{signature(con = "CNSet")}: ... }
    \item{coerce}{\code{signature(from="CNSetLM")}: ... }
    \item{corr}{\code{signature(object = "CNSet", allele = "character")}: ... }
    \item{flags}{\code{signature(object="CNSet")}: SNP flags }
    \item{initialize}{\code{signature(.Object = "CNSet")}: ... }
%    \item{lM}{\code{signature(object = "CNSet")}: ... }
%    \item{lM<-}{\code{signature(object = "CNSet", value = "LinearModelParameter")}: ... }
    \item{nu}{\code{signature(object = "CNSet", allele = "character")}: ... }
    \item{open}{\code{signature(con = "CNSet")}: ... }
    \item{phi}{\code{signature(object = "CNSet", allele = "character")}: ... }
    \item{sigma2}{\code{signature(object = "CNSet", allele = "character")}: ... }
    \item{tau2}{\code{signature(object = "CNSet", allele = "character")}: ... }
	 }
}

\author{
R. Scharpf
}

\seealso{
}
\examples{
 	if(require("genomewidesnp6Crlmm")){
		require("genomewidesnp6Crlmm")
		fns <- c("SNP_A-2131660", "SNP_A-1967418", "SNP_A-1969580", "SNP_A-4263484",
			 "SNP_A-1978185", "SNP_A-4264431", "SNP_A-1980898", "SNP_A-1983139",
			 "SNP_A-4265735", "SNP_A-1995832")
		theCalls <- matrix(2, nc=2, nrow=10)
		A <- matrix(sample(1:1000, 20), 10,2)
		B <- matrix(sample(1:1000, 20), 10,2)
		p <- matrix(runif(20), nc=2)
		theConfs <- round(-1000*log2(1-p))
		## Batch can be defined by the scan date of the array
		##or the 96 well chemistry plate from which the
		##samples were derived.  Here we indicate that the two
		##samples were from the same batch.
		batch <- rep(factor(1), ncol(A))
		## each parameter is a R x C matrix, where the number
		## of rows (R) corresponds to the number of features
		## and the number of columns (C) corresponds to the
		## number of batches.  In this toy example, the
		## samples were assumed to be from the same batch.
		## Ordinarily, one would have 50+ samples in a given
		## batch.
		dns <- list(fns, batch="1")
		obj <- new("CNSet",
			   alleleA=A,
			   alleleB=B,
			   call=theCalls,
			   callProbability=theConfs,
			   batch=as.factor(rep(1, ncol(A))),
			   annotation="genomewidesnp6")
		assayDataElementNames(batchStatistics(obj))
		featureNames(obj) <- fns
		## Accessors
		calls(obj)
		confs(obj)
		A(obj)
		B(obj)
		featureData(obj) <- addFeatureAnnotation(obj)
		isSnp(obj)
		chromosome(obj)
		position(obj)
	}

}
\keyword{classes}
