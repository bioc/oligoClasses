setValidity("CNSet", function(object) assayDataValidMembers(assayData(object), c("alleleA", "alleleB", "call", "callProbability")))
setMethod("A", "CNSet", function(object, ...) allele(object, "A", ...))
setMethod("B", "CNSet", function(object, ...) allele(object, "B", ...))
setReplaceMethod("A", "CNSet", function(object, value) {
	assayDataElementReplace(object, "CNSet", value)
})
setReplaceMethod("B", "CNSet", function(object, value) {
	assayDataElementReplace(object, "CNSet", value)
})
setMethod("close", "CNSet", function(con, ...){
	##con is just to keep the same generic arguments
	object <- con
	if(!isFF(object)) return()
	names <- ls(assayData(object))
	L <- length(names)
	for(i in 1:L) close(eval(substitute(assayData(object)[[NAME]], list(NAME=names[i]))))
	return()
})

setMethod("open", "CNSet", function(con, ...){
	object <- con
	if(!isFF(object)) return()
	names <- ls(assayData(object))
	L <- length(names)
	for(i in 1:L) open(eval(substitute(assayData(object)[[NAME]], list(NAME=names[i]))))
	return()
})
	    





