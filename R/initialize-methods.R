setMethod("initialize", "oligoSnpSet",
	  function(.Object,
		   assayData=assayDataNew(call=call,
		                          callProbability=callProbability,
		                          copyNumber=copyNumber,
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
		  browser()
		  if(missing(assayData)){
			  .Object <- callNextMethod(.Object,
						    call=call,
						    callProbability=callProbability,
						    copyNumber=copyNumber,
						    cnConfidence=cnConfidence,...)
		  } else{
			  .Object <- callNextMethod(.Object,
						    assayData=assayData, ...)
		  }
		  if(missing(annotation)){
			  stop("must specify annotation")
		  } else{
			  if(!isSupportedAnnotation(annotation) & (missing(position) | missing(chromosome) | missing(isSnp))){
				  stop("annotation is not supported. Must specify chromosome, position, and an indicator for whether the marker is polymorphic (isSnp)")
			  } else{
				  pData(featureData)$chromosome <- chromosome
				  pData(featureData)$position <- position
				  pData(featureData)$isSnp <- isSnp
			  }
			  .Object@featureData <- featureData
		  }
		  ## Do after annotation has been assigned
		  if(!(all(c("chromosome", "position", "isSnp") %in% varLabels(featureData))) & isSupportedAnnotation(annotation)){
			  .Object@featureData <- addFeatureAnnotation.SnpSet(.Object)
		  }
		  return(.Object)
	  })

isSupportedAnnotation <- function(x){
	validAnn <- validAnnotation()
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
