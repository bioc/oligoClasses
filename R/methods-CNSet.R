setMethod("CA", "CNSet", function(object) assayData(object)[["CA"]]/100)
setMethod("CB", "CNSet", function(object) assayData(object)[["CB"]]/100)
setReplaceMethod("CA", signature(object="CNSet", value="matrix"),
		 function(object, value){
			 value <- matrix(as.integer(value*100), nrow(value), ncol(value), dimnames=dimnames(value))
			 assayDataElementReplace(object, "CA", value)		
		 })
setReplaceMethod("CB", signature(object="CNSet", value="matrix"),
		 function(object, value){
			 value <- matrix(as.integer(value*100), nrow(value), ncol(value), dimnames=dimnames(value))			 
			 assayDataElementReplace(object, "CB", value)
		 })

##setMethod("rangedData", "CNSet", function(object) segmentData(object))
##setReplaceMethod("rangedData", c("CNSet", "RangedData"), function(object, value){
##	segmentData(object) <- value
##})
##
##setMethod("segmentData", "CNSet", function(object) object@segmentData)
##setReplaceMethod("segmentData", c("CNSet", "RangedData"), function(object, value){
##	object@segmentData <- value
##	object
##})
##
##setMethod("emissionPr", "CNSet", function(object) object@emissionPr)
##setReplaceMethod("emissionPr", c("CNSet", "array"), function(object, value){
##	object@emissionPr <- value
##	object
##})

setMethod("show", "CNSet", function(object){
	callNextMethod(object)
##	cat("emissionPr\n")
##	cat("   array:", nrow(object), "features,", ncol(object), "samples,", dim(emissionPr(object))[3], "states\n")
##	cat("   hidden states:\n")
##	cat("      ", dimnames(emissionPr(object))[[3]], "\n")
##	cat("   Missing values:", sum(is.na(emissionPr(object))), "\n")
##	if(!all(is.na(emissionPr(object)))){
##		cat("   minimum value:", min(emissionPr(object), na.rm=TRUE), "\n")
##	} else  cat("   minimum value: NA (all missing)\n")
##	cat("rangedData:  ")
##	cat("    ", show(rangedData(object)), "\n")
})

##setMethod("start", "CNSet", function(x, ...) start(segmentData(x), ...))
##setMethod("end", "CNSet", function(x, ...) end(segmentData(x), ...))
##setMethod("width", "CNSet", function(x) width(segmentData(x)))
