setMethod("show", "CNSet", function(object){
	callNextMethod(object)
	cat("lM: ", length(lM(object)), " parameters, ", nrow(object), " features, ", length(unique(batch(object))), " batches\n")
	cat("   element names: ", paste(ls(lM(object))[1:4], collapse=",  "), "\n")
	cat("                  ", paste(ls(lM(object))[5:8], collapse=",  "), "\n")
	cat("                  ", paste(ls(lM(object))[9:12],collapse=",  "), "\n")
	cat("                  ", paste(ls(lM(object))[13],  collapse=",  "),   "\n")
})


setMethod("[", "CNSet", function(x, i, j, ..., drop=FALSE){
	x <- callNextMethod(x, i, j, ..., drop=drop)
	if(!missing(i)){
		if(class(lM(x)) == "ffdf"){
			lM(x) <- lapply(physical(lM(x)), function(x, i){open(x); x[i, ]}, i=i)
		} else {
			lM(x) <- lapply(lM(x), function(x, i) x[i, , drop=FALSE], i=i)
		}
	}
	x
})
setMethod("batch", "CNSet", function(object) object@batch)

setMethod("batchNames", "CNSet", function(object)  batchNames(lM(object)))

setReplaceMethod("batchNames", "CNSet", function(object, value) {
	batchNames(object@lM) <- value
	return(object)
})

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
setReplaceMethod("lM", signature=signature(object="CNSet", value="LinearModelParameter"),
		 function(object, value){
			 object@lM <- value
			 object
		 })

