setMethod("show", "CNSet", function(object){
	callNextMethod(object)
	bns <- batchNames(object)
	freq <- as.integer(table(batch(object)))
	cat("batch:   ", paste(bns, freq, collapse=", "), "\n")
	cat("lM: ", length(lM(object)), " parameters, ", nrow(object), " features, ", length(unique(batch(object))), " batches\n")
	cat("   element names: ", paste(ls(lM(object))[1:4], collapse=",  "), "\n")
	cat("                  ", paste(ls(lM(object))[5:8], collapse=",  "), "\n")
	cat("                  ", paste(ls(lM(object))[9:12],collapse=",  "), "\n")
	cat("                  ", paste(ls(lM(object))[13],  collapse=",  "),   "\n")
})

setMethod("[", "CNSet", function(x, i, j, ..., drop=FALSE){
	x <- callNextMethod(x, i, j, ..., drop=drop)
	if(!missing(j)){
		batch(x) <- batch(x)[j]
		nms <- unique(as.character(batch(x)))
		## need to subset columns of LinearModelParameter
		## Adapted from the '[' method for eSet in Biobase
		## redefine 'j'
		j <- nms %in% batchNames(x)
		storage.mode <- storageMode(lM(x))
		## i (if defined) is already subset by callNextMethod
		orig <- lM(x)
		lM(x) <-
			switch(storage.mode,
			       environment =,
			       lockedEnvironment = {
				       aData <- new.env(parent=emptyenv())
				       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i, j, ..., drop = drop]
				       if ("lockedEnvironment" == storage.mode) Biobase:::assayDataEnvLock(aData)
				       aData
			       },
			       list = {
				       lapply(orig, function(obj) obj[i, j, ..., drop = drop])
			       })
	}
	x
})
setMethod("batch", "CNSet", function(object) object@batch)

setReplaceMethod("batch", signature=signature(object="CNSet"),
	 function(object, value){
		 object@batch <- as.factor(value)
		 object
})

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

setMethod("nu", c("CNSet", "character"), function(object, allele) nu(lM(object), allele))
setMethod("phi", c("CNSet", "character"), function(object, allele) phi(lM(object), allele))
setMethod("sigma2", c("CNSet", "character"), function(object, allele) sigma2(lM(object), allele))
setMethod("tau2", c("CNSet", "character"), function(object, allele) tau2(lM(object), allele))
setMethod("corr", c("CNSet", "character"), function(object, allele) corr(lM(object), allele))

