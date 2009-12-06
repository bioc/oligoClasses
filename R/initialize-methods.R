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

setValidity("oligoSnpSet", function(object) {
	assayDataValidMembers(assayData(object), c("call", "callProbability", "copyNumber", "cnConfidence"))
})

## RS: ask BC about this... initialization method for CNSet does not work when this is uncommented
##setValidity("AlleleSet",
##            function(object){
##              grp1 <- c("alleleA", "alleleB")
##              grp2 <- c("senseAlleleA", "senseAlleleB",
##                        "antisenseAlleleA", "antisenseAlleleB")
##              elem <- assayDataElementNames(object)
##              ok <- all(grp1 %in% elem) || all(grp2 %in% elem)
##              f <- function(x) paste("'", x, "'", collapse=" + ", sep="")
##              if (!ok){
##                paste("Elements of 'AlleleSummarySet' must be:",
##                      f(grp1), "OR", f(grp2))
##              }else{
##                TRUE
##              }
##            })

setMethod("initialize", "CNSet",
          function(.Object,
		   call=new("matrix"),		   		   
		   CA=matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
		   CB=matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
                   callProbability=matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
                   alleleA = matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
                   alleleB = matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
                   phenoData = annotatedDataFrameFrom(call, byrow=FALSE),		   
		   featureData=annotatedDataFrameFrom(call, byrow=TRUE),
		   experimentData=new("MIAME"),
		   protocolData=phenoData[, integer(0)],
		   segmentData=new("RangedData"),
		   emissionPr=new("array"),
		   position=integer(),
		   chromosome=integer(),
		   isSnp=integer(),
		   annotation=character(), ... ){
		  ##The ... should be additional assayDataElements (e.g., for a class that extends CNSet)
		  .Object <- callNextMethod(.Object,
					    call=call,
					    callProbability=callProbability,
					    alleleA=alleleA,
					    alleleB=alleleB,
					    CA=CA,
					    CB=CB,
					    phenoData=phenoData,
					    featureData=featureData,
					    experimentData=experimentData,
					    protocolData=protocolData,
					    annotation=annotation, ...)
		  annotation <- .Object@annotation
		  ##add chromosome, position, isSnp to featureData
		  message("Adding chromosome, position, isSnp indicator to featureData...")
		  if(length(annotation) < 1){
			  if((length(position) < 1| length(chromosome) < 1 | length(isSnp) < 1)){
				  stop("must specify annotation if 'chromosome', 'position', and 'isSnp' are missing")
			  } else {
				  pData(featureData)$chromosome <- chromosome
				  pData(featureData)$position <- position
				  pData(featureData)$isSnp <- isSnp
			  }
		  } else{
			  if((length(position) < 1| length(chromosome) < 1 | length(isSnp) < 1)){
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
		  ##Do after annotation has been assigned
		  if(!(all(c("chromosome", "position", "isSnp") %in% varLabels(featureData))) & isSupportedAnnotation(annotation)){
			  .Object@featureData <- addFeatureAnnotation(.Object)
		  }		  
		  if(!missing(emissionPr)) .Object@emissionPr <- emissionPr
		  segmentData(.Object) <- segmentData
		  .Object	    
          })


setValidity("CNSet", function(object) {
	assayDataValidMembers(assayData(object), c("CA", "CB", "call", "callProbability", "alleleA", "alleleB"))
})
