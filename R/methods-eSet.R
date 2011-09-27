setMethod("annotate", "eSet", function(object){
	annotation <- object@annotation
	featureData <- object@featureData
	position <- grep("position", varLabels(featureData))
	chromosome <- grep("chromosome", varLabels(featureData))
	isSnp <- grep("isSnp", varLabels(featureData))
	if(length(annotation) < 1){
		if((length(position) < 1| length(chromosome) <1 | length(isSnp) <1)){
			stop("must specify annotation if 'chromosome', 'position', and 'isSnp' are missing")
		} else {
			pData(featureData)$chromosome <- chromosome
			pData(featureData)$position <- position
			pData(featureData)$isSnp <- isSnp
		}
	} else{
		if((length(position) < 1| length(chromosome) < 1| length(isSnp) < 1)){
			if(!isSupportedAnnotation(annotation)){
				stop("The annotation is not supported. Arguments 'chromosome', 'position', and 'isSnp' can be omitted from the initialization only if the annotation is supported (see oligoClasses:::supportedAnnotation()).")
			}
		} else {
			pData(featureData)$chromosome <- chromosome
			pData(featureData)$position <- position
			pData(featureData)$isSnp <- isSnp
		}
		object@featureData <- featureData
	}
	## Do after annotation has been assigned
	if(!(all(c("chromosome", "position", "isSnp") %in% varLabels(featureData))) & isSupportedAnnotation(annotation)){
		object@featureData <- addFeatureAnnotation(object)
	}
	return(object)
})

setMethod("isSnp", "eSet", function(object) {
	labels <- fvarLabels(object)
	if("isSnp" %in% labels){
		res <- fData(object)[, "isSnp"]
	} else{
		res <- as.integer(featureNames(object) %in% snpNames(object))
	}
	return(res==1)
})


setMethod("db", "eSet",
          function(object) {
		  requireAnnotation(annotation(object)) || stop(paste(annotation(object), "package not available"))
		  get(annotation(object))@getdb()
	  })

setMethod("chromosome", "eSet",
	  function(object, na.rm=FALSE){
		  if(!("chromosome" %in% fvarLabels(object))){
			  stop("chromosome not in fvarLabels")
		  }
		  chrom <- chromosome(featureData(object), na.rm)
		  return(chrom)
	  })

setReplaceMethod("chromosome", "eSet",
		 function(object, value){
			 fData(object)$chromosome <-  value
			 object
		 })


setMethod("position", "eSet",
          function(object, na.rm=FALSE){
		  if(!("position" %in% fvarLabels(object))){
			  stop("position not in fvarLabels")
		  }
		  pos <- position(featureData(object), na.rm)
		  return(pos)
          })
