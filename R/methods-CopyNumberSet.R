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



##setMethod("order", "CopyNumberSet",
##	  function(object, ..., na.last=TRUE, decreasing=FALSE){
##		  chromosomePositionOrder(...)
##	  })
