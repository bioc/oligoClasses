setMethod("batchNames", "LinearModelParameter",
	  function(object){
		  ##should call method for AssayData
		  sampleNames(object)
	  })

setReplaceMethod("batchNames", "LinearModelParameter",
		 function(object, value){
			 sampleNames(object) <- value
			 return(object)
		 })

