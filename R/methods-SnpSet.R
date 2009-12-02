setMethod("calls", "SnpSet", function(object) assayData(object)$call)
setReplaceMethod("calls", signature(object="SnpSet", value="matrix"), function(object, value) assayDataElementReplace(object, "call", value))
setMethod("confs", "SnpSet", function(object) {
	X <- assayData(object)$callProbability
	P <- 1-exp(-X/1000)
	return(P)
})
setReplaceMethod("confs", signature(object="SnpSet", value="matrix"),
		 function(object, value){
			 ##convert probability to integer
			 P <- value
			 dns <- dimnames(P)
			 X <- -1000*log(1-P)
			 X <- matrix(as.integer(X), nrow(X), ncol(X))
			 dimnames(X) <- dns
			 assayDataElementReplace(object, "callProbability", X)
		 })

setMethod("callsConfidence", "SnpSet", function(object)
          assayDataElement(object, "callProbability"))

setReplaceMethod("callsConfidence", signature(object="SnpSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "callProbability", value))
