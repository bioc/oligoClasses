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
	  function(object){
		  if(!("chromosome" %in% fvarLabels(object))){
			  stop("chromosome not in fvarLabels")
		  } 
		  return(featureData(object)$chromosome)
	  })

setReplaceMethod("chromosome", c("eSet", "ANY"),
		 function(object, value){
			 fData(object)$chromosome <-  value
			 object
		 })


setMethod("position", "eSet",
          function(object){
		  if(!("position" %in% fvarLabels(object))){
			  stop("position not in fvarLabels")
		  }
		  pos <- featureData(object)$position
          })
