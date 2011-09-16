setMethod("copyNumber", "oligoSnpSet", function(object) {
	cn <- assayDataElement(object, "copyNumber")
	if(vmode(cn) == "double") {
		return(cn)
	}
	if(vmode(cn) == "integer"){
		if(is(cn, "ff")) return(cn)
		cn <- cn/100
		return(cn)
	}
})


setReplaceMethod("copyNumber", signature(object="oligoSnpSet", value="matrix"),
                 function(object, value){
			 ##value <- matrix(as.integer(value*100), nrow(value), ncol(value), dimnames=dimnames(value))
			 assayDataElementReplace(object, "copyNumber", value)
		 })

setMethod("cnConfidence", "oligoSnpSet", function(object) assayData(object)[["cnConfidence"]])
setReplaceMethod("cnConfidence", signature(object="oligoSnpSet", value="matrix"),
                 function(object, value){
##			 conf <- value
##			 dns <- dimnames(conf)
##			 X <- matrix(as.integer(conf*100), nrow(conf), ncol(conf))
##			 dimnames(X) <- dns
			 assayDataElementReplace(object, "cnConfidence", value)
                 })
