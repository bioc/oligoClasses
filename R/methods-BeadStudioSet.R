setMethod("updateObject", signature(object="BeadStudioSet"),
          function(object, ..., verbose=FALSE) {
		  if (verbose) message("updateObject(object = 'BeadStudioSet')")
		  obj <- tryCatch(callNextMethod(object), error=function(e) NULL)
		  if(is.null(obj)){
			  obj <- new("BeadStudioSet",
				     assayData = updateObject(assayData(object),
				     ..., verbose=verbose),
				     phenoData = phenoData(object),
				     experimentData = updateObject(experimentData(object),
				     ..., verbose=verbose),
				     annotation = updateObject(annotation(object),
				     ..., verbose=verbose),
				     featureData=updateObject(featureData(object), ..., verbose=FALSE),
				     ...)
		  }
		  if (all(isCurrent(obj))) return(obj)
		  obj
          })

setMethod("lrr", "BeadStudioSet", function(object){
	return(assayDataElement(object, "lrr"))
})

setMethod("copyNumber", "BeadStudioSet", function(object)
	  lrr(object))

setReplaceMethod("lrr", c("BeadStudioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "lrr", value)
	 })

setReplaceMethod("lrr", c("BafLrrSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "lrr", value)
	 })

setReplaceMethod("copyNumber", c("BeadStudioSet", "ANY"),
		 function(object, value) {
			 lrr(object) <- value
			 object
	 })

setMethod("baf", "BeadStudioSet",
	  function(object) {
		  return(assayDataElement(object, "baf"))
	  })

setReplaceMethod("baf", c("BeadStudioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "BAF", value)
	 })

setAs("BeadStudioSet", "data.frame",
      function(from, to){
	      cn <- as.numeric(lrr(from))/100
	      bf <- as.numeric(baf(from))/1000
	      x <- rep(position(from)/1e6, ncol(from))
	      ##x <- rep(position(object)[marker.index], 4)/1e6
	      is.snp <- rep(isSnp(from), ncol(from))
	      id <- rep(sampleNames(from), each=nrow(from))
	      df <- data.frame(x=x, lrr=cn, baf=bf, id=id,
			       is.snp=is.snp,
			       stringsAsFactors=FALSE)
	      df$id <- factor(df$id, ordered=TRUE, levels=unique(df$id))
	      return(df)
      })

setMethod("show", signature(object="BeadStudioSet"), function(object){
	callNextMethod(object)
	##cat("Genome Build: ", genomeBuild(object), "\n")
	##cat("Integer representation of BAF/LRR: ", isInteger(object), "\n")
})
