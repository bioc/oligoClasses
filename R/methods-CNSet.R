##setMethod("CA", "CNSet", function(object) assayData(object)[["CA"]]/100)
##setMethod("CB", "CNSet", function(object) assayData(object)[["CB"]]/100)
setMethod("CB", "CNSet", function(object) assayDataElement(object, "CB"))
setMethod("CA", "CNSet", function(object) assayDataElement(object, "CA"))
setReplaceMethod("CB", "CNSet", function(object, value) assayDataElementReplace(object, "CB", value))
setReplaceMethod("CA", "CNSet", function(object, value) assayDataElementReplace(object, "CA", value))
##setReplaceMethod("CA", signature(object="CNSet", value="matrix"),
##		 function(object, value){
##			 value <- matrix(as.integer(value*100), nrow(value), ncol(value), dimnames=dimnames(value))
##			 assayDataElementReplace(object, "CA", value)		
##		 })
##setReplaceMethod("CB", signature(object="CNSet", value="matrix"),
##		 function(object, value){
##			 value <- matrix(as.integer(value*100), nrow(value), ncol(value), dimnames=dimnames(value))			 
##			 assayDataElementReplace(object, "CB", value)
##		 })


