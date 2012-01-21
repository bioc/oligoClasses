setMethod("lrr", "BeadStudioSet", function(object)
	  assayDataElement(object, "lrr"))
setReplaceMethod("lrr", c("BeadStudioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "lrr", value)
	 })
setMethod("baf", "BeadStudioSet",
	  function(object) {
		  assayDataElement(object, "baf")
	 })
setReplaceMethod("baf", c("BeadStudioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "BAF", value)
	 })
setAs("BeadStudioSet", "data.frame",
      function(from, to){
	      cn <- as.numeric(lrr(from))
	      bf <- as.numeric(baf(from))
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
