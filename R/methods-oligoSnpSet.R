setMethod("copyNumber", "oligoSnpSet", function(object){
	X <- assayData(object)[["copyNumber"]]
	CN <- X/100
	return(CN)
})
setReplaceMethod("copyNumber", signature(object="oligoSnpSet", value="matrix"),
                 function(object, value){
			 CN <- value
			 dns <- dimnames(CN)
			 X <- matrix(as.integer(CN*100), nrow(CN), ncol(CN))
			 dimnames(X) <- dns
			 assayDataElementReplace(object, "copyNumber", X)
		 })
setMethod("cnConfidence", "oligoSnpSet", function(object){
	X <- assayData(object)[["cnConfidence"]]
	conf <- X/100
	return(X)	
  })
setReplaceMethod("cnConfidence", signature(object="oligoSnpSet", value="matrix"),
                 function(object, value){
			 conf <- value
			 dns <- dimnames(conf)
			 X <- matrix(as.integer(conf*100), nrow(conf), ncol(conf))
			 dimnames(X) <- dns
			 assayDataElementReplace(object, "cnConfidence", X)
                 })
