setValidity("CNSet", function(object) assayDataValidMembers(assayData(object), c("alleleA", "alleleB", "call", "callProbability")))
setMethod("allele", "CNSet",
          function(object, allele){
            stopifnot(!missing(allele))
            allele <- match.arg(allele, c("A", "B"))
	    what <- paste("allele", allele, sep="")
            assayDataElement(object, what)
          })
setMethod("A", "CNSet", function(object, ...) allele(object, "A", ...))
setMethod("B", "CNSet", function(object, ...) allele(object, "B", ...))
setReplaceMethod("A", "CNSet", function(object, value) {
	assayDataElementReplace(object, "alleleA", value)
})
setReplaceMethod("B", "CNSet", function(object, value) {
	assayDataElementReplace(object, "alleleB", value)
})
setMethod("close", "CNSet", function(con, ...){
	##con is just to keep the same generic arguments
	object <- con
	if(!isFF(object)) return()
	names <- ls(assayData(object))
	L <- length(names)
	for(i in 1:L) close(eval(substitute(assayData(object)[[NAME]], list(NAME=names[i]))))
	physical <- get("physical")
	lapply(physical(lM(con)), open)
	return()
})

setMethod("open", "CNSet", function(con, ...){
	object <- con
	if(!isFF(object)) return()
	names <- ls(assayData(object))
	L <- length(names)
	for(i in 1:L) open(eval(substitute(assayData(object)[[NAME]], list(NAME=names[i]))))
	physical <- get("physical")
	lapply(physical(lM(con)), close)	
	return()
})
	    
setMethod("lM", "CNSet", function(object) object@lM)
setReplaceMethod("lM", c("CNSet", "list_or_ffdf"), function(object, value){
	object@lM <- value
	object
})




