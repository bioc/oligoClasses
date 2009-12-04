setMethod("initialize", "oligoSnpSet",
	  function(.Object,
		   assayData=assayDataNew(call=call,
		                          callProbability=matrix(as.integer(-1000*(log(1-callProbability))), nrow(call), ncol(call), dimnames=dimnames(call)),
		                          copyNumber=matrix(as.integer(copyNumber*100), nrow(call), ncol(call), dimnames=dimnames(call)),
                                          cnConfidence=cnConfidence, ...),
		   call=new("matrix"),
		   callProbability=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
		   copyNumber=matrix(numeric(), nrow=nrow(call), ncol=ncol(call),  dimnames=dimnames(call)),
		   cnConfidence=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
		   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
		   position,
		   chromosome,
		   isSnp,
		   annotation, ... ){
		  message("Storing copyNumber and callProbability as integers.  Use the copyNumber() and confs() accessors")
		  .Object <- callNextMethod(.Object,
					    assayData=assayData,
					    featureData=featureData, ...)
		  if(missing(annotation)){
			  if((missing(position) | missing(chromosome) | missing(isSnp))){
				  stop("must specify annotation if 'chromosome', 'position', and 'isSnp' are missing")
			  } else {
				  pData(featureData)$chromosome <- chromosome
				  pData(featureData)$position <- position
				  pData(featureData)$isSnp <- isSnp
			  }
		  } else{
			  .Object@annotation <- annotation
			  if((missing(position) | missing(chromosome) | missing(isSnp))){
				  if(!isSupportedAnnotation(annotation)){
					  stop("The annotation is not supported. Arguments 'chromosome', 'position', and 'isSnp' can be omitted from the initialization only if the annotation is supported (see oligoClasses:::supportedAnnotation()).")
				  }
			  } else {
				  pData(featureData)$chromosome <- chromosome
				  pData(featureData)$position <- position
				  pData(featureData)$isSnp <- isSnp
			  }
			  .Object@featureData <- featureData
		  }
		  ## Do after annotation has been assigned
		  if(!(all(c("chromosome", "position", "isSnp") %in% varLabels(featureData))) & isSupportedAnnotation(annotation)){
			  .Object@featureData <- addFeatureAnnotation(.Object)
		  }
		  return(.Object)
	  })

isSupportedAnnotation <- function(x){
	validAnn <- supportedAnnotation()
	x %in% validAnn
}

supportedAnnotation <- function(){
	c("pd.mapping50k.hind240", "pd.mapping50k.xba240",
	  "pd.mapping50k.hind240,pd.mapping50k.xba240",	  
	  "pd.mapping250k.nsp",
	  "pd.mapping250k.sty",
	  "pd.mapping250k.nsp,pd.mapping250k.sty",
	  "pd.genomewidesnp.5",
	  "pd.genomewidesnp.6",
	  "genomewidesnp6",
	  "genomewidesnp5",
	  "human370v1c",
	  "human370quadv3c",
	  "human550v3b",
	  "human650v3a",
	  "human610quadv1b",
	  "human660quadv1a",
	  "human1mduov3b")	  
}

setValidity("oligoSnpSet", function(object) {
	assayDataValidMembers(assayData(object), c("call", "callProbability", "copyNumber", "cnConfidence"))
})

setValidity("AlleleSet",
            function(object){
              grp1 <- c("alleleA", "alleleB")
              grp2 <- c("senseAlleleA", "senseAlleleB",
                        "antisenseAlleleA", "antisenseAlleleB")
              elem <- assayDataElementNames(object)
              ok <- all(grp1 %in% elem) || all(grp2 %in% elem)
              f <- function(x) paste("'", x, "'", collapse=" + ", sep="")
              if (!ok){
                paste("Elements of 'AlleleSummarySet' must be:",
                      f(grp1), "OR", f(grp2))
              }else{
                TRUE
              }
            })
