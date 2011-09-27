##setMethod("copyNumber", "CopyNumberSet", function(object) assayData(object)[["copyNumber"]])
setMethod("copyNumber", "CopyNumberSet", function(object) {
	cn <- assayDataElement(object, "copyNumber")
	return(cn)
})

setReplaceMethod("copyNumber", signature(object="CopyNumberSet", value="matrix"),
                 function(object, value){
			 assayDataElementReplace(object, "copyNumber", value)
		 })

setMethod("cnConfidence", "CopyNumberSet", function(object) assayData(object)[["cnConfidence"]])
setReplaceMethod("cnConfidence", signature(object="CopyNumberSet", value="matrix"),
                 function(object, value){
			 assayDataElementReplace(object, "cnConfidence", value)
                 })

setMethod("checkOrder", signature(object="CopyNumberSet"),
	  function(object, verbose=FALSE){
		  .checkOrder(object, verbose)
	  })

setMethod("order", "CopyNumberSet",
	  function(..., na.last=TRUE, decreasing=FALSE){
		  object <- list(...)[[1]]
		  chromosomePositionOrder(object, ...)
	  })
