\name{CNSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CNSet-class}
\alias{initialize,CNSet-method}
\alias{allele,CNSet-method}
\alias{A,CNSet-method}
\alias{B,CNSet-method}
\alias{A<-,CNSet-method}
\alias{B<-,CNSet-method}
\alias{close,CNSet-method}
\alias{open,CNSet-method}
	       
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
		 \code{genotype} function.

		 The following creates a very basic \code{CNSet} with
		 \code{assayData} containing the required alleleA,
		 alleleB, call, and callProbability matrices.

		 \code{new(CNSet, alleleA=new("matrix"), alleleB=new("matrix"), call=new("matrix"),
			   callProbability=new("matrix"))}
		 
}
\section{Slots}{
  \describe{
	  \item{\code{lM}:}{Object of class \code{"list_or_ffdf"}. For storing parameters related to copy number estimation.}
	  \item{\code{assayData}:}{Object of class \code{"AssayData"} }
	  \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} }
	  \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} }
	  \item{\code{experimentData}:}{Object of class \code{"MIAME"} }
	  \item{\code{annotation}:}{Object of class \code{"character"}}
	  \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} }
	  \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
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
  \item{A}{\code{signature(object = "CNSet")}: accessor for normalized intensities for A allele}
  \item{B}{\code{signature(object = "CNSet")}: accessor for normalized intensities for B allele }
  \item{calls}{\code{signature(object = "CNSet")}: accessor for crlmm genotype calls}
  \item{confs}{\code{signature(object = "CNSet")}: accessor for genotype confidence scores}
}
}
\author{R. Scharpf}

\seealso{
	\code{\linkS4class{AlleleSet}}, \code{\linkS4class{SnpSuperSet}}
}
\examples{
	##An empty CNSet
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
		##lM <- crlmm:::initializeParamObject(dns)
		obj <- new("CNSet",
			   alleleA=A,
			   alleleB=B,
			   call=theCalls,
			   callProbability=theConfs,
			   lM=list(),
			   annotation="genomewidesnp6")
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
