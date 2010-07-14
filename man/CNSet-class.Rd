\name{CNSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CNSet-class}
\alias{CA}
\alias{CA<-}
\alias{CA,CNSet-method}
\alias{CA<-,CNSet,matrix-method}
\alias{CA<-,CNSet-method}
\alias{CB}
\alias{CB<-}
\alias{CB,CNSet-method}
\alias{CB<-,CNSet-method}
\alias{CB<-,CNSet,matrix-method}
\alias{initialize,CNSet-method}
\alias{totalCopyNumber,CNSet,integer,integer-method}

\title{Class "CNSet"}
\description{A class to store copy number estimates as well as locus-level summaries of the normalized intensities}
\section{Objects from the Class}{
		 Objects can be created by calls of the form 
		 \code{new("CNSet", call, CA, CB, callProbability, alleleA, alleleB, phenoData, featureData, experimentData, protocolData, segmentData, emissionPr, position, chromosome, isSnp, annotation, ...)}.

}
\section{Slots}{
  \describe{
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
Class \code{"\linkS4class{SnpSuperSet}"}, directly.
Class \code{"\linkS4class{AlleleSet}"}, by class "SnpSuperSet", distance 2.
Class \code{"\linkS4class{SnpSet}"}, by class "SnpSuperSet", distance 2.
Class \code{"\linkS4class{eSet}"}, by class "SnpSuperSet", distance 3.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "SnpSuperSet", distance 4.
Class \code{"\linkS4class{Versioned}"}, by class "SnpSuperSet", distance 5.
}
\section{Methods}{
  \describe{
    \item{CA}{\code{signature(object = "CNSet")}: ... }
    \item{CA<-}{\code{signature(object = "CNSet", value = "matrix")}: ... }
    \item{CB}{\code{signature(object = "CNSet")}: ... }
    \item{CB<-}{\code{signature(object = "CNSet", value = "matrix")}: ... }
   \item{initialize}{\code{signature(.Object = "CNSet")}: ... }
  \item{totalCopyNumber}{\code{signature(object="CNSet", i="integer", j="integer")}: Computes total copy number for markers \code{i} and samples \code{j} of a \code{CNSet} object}
}
}
\author{R. Scharpf}

\seealso{
	\code{\linkS4class{AlleleSet}}, \code{\linkS4class{SnpSuperSet}}
}
\examples{
	##For purposes of illustration only
	theCalls <- matrix(2, nc=2, nrow=10)
	A <- matrix(sample(1:1000, 20), 10,2)
	B <- matrix(sample(1:1000, 20), 10,2)
        ## stored CA, CB as integers:
	CA <- matrix(as.integer(rnorm(20, 100)), nrow=10)
	CB <- matrix(as.integer(rnorm(20, 100)), nrow=10)
	p <- matrix(runif(20), nc=2)
	theConfs <- round(-1000*log2(1-p))
	obj <- new("CNSet",
		   alleleA=A,
		   alleleB=B,
		   call=theCalls,
		   callProbability=theConfs,
		   CA=CA,
		   CB=CB)

	## Accessors
	calls(obj)
	confs(obj)
	A(obj)
	B(obj)
	CA(obj)
	CB(obj)

	##with annotation
	if(require("genomewidesnp6Crlmm")){
		ids <- c("SNP_A-2131660", "SNP_A-1967418", "SNP_A-1969580", "SNP_A-4263484",
			 "SNP_A-1978185", "SNP_A-4264431", "SNP_A-1980898", "SNP_A-1983139",
			 "SNP_A-4265735", "SNP_A-1995832")
		rownames(theCalls) <- rownames(p) <- rownames(theConfs) <- ids
		rownames(A) <- rownames(B) <- rownames(CA) <- rownames(CB) <- ids
		obj2 <- new("CNSet",
			   alleleA=A,
			   alleleB=B,
			   call=theCalls,
			   callProbability=theConfs,
			   CA=CA,
			   CB=CB, annotation="genomewidesnp6")
		featureData(obj2) <- addFeatureAnnotation(obj2)	   
		fvarLabels(obj2)
		## additional accessors
		isSnp(obj2)
		chromosome(obj2)
		position(obj2)
		total.copy.number <- totalCopyNumber(obj2, i=1:5, j=1:2)
		}
}
	
\keyword{classes}
