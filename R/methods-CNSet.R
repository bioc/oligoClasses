setMethod("show", "CNSet", function(object){
	callNextMethod(object)
	bns <- batchNames(object)
	freq <- as.integer(table(batch(object)))
	cat("batch:   ", paste(bns, freq, collapse=", "), "\n")
	cat("lM: ", length(lM(object)), " parameters, ", nrow(object), " features, ", length(unique(batch(object))), " batches\n")
	cat("   element names: ", paste(ls(lM(object))[1:4], collapse=",  "), "\n")
	cat("                  ", paste(ls(lM(object))[5:8], collapse=",  "), "\n")
	cat("                  ", paste(ls(lM(object))[9:12],collapse=",  "), "\n")
	cat("                  ", paste(ls(lM(object))[13:length(lM(object))],  collapse=",  "),   "\n")
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

##assayDataElement <- function(object, elt) assayData(object)[[elt]]
##assayDataElementReplace <- function(obj, elt, value) {
##    storage.mode <- storageMode(obj)
##    switch(storageMode(obj),
##           "lockedEnvironment" = {
##               aData <- copyEnv(assayData(obj))
##               if (is.null(value)) rm(list=elt, envir=aData)
##               else aData[[elt]] <- value
##               assayDataEnvLock(aData)
##               assayData(obj) <- aData
##           },
##           "environment" = {
##               if (is.null(value)) rm(list=elt, envir=assayData(obj))
##               else assayData(obj)[[elt]] <- value
##           },
##           list = assayData(obj)[[elt]] <- value)
##    obj
##}

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
	open(assayDataElement(object, "alleleA")) ## should do nothing if matrix
	obj <- assayDataElementReplace(object, "alleleA", value)
	close(assayDataElement(object, "alleleA")) ## should do nothing if matrix
})

setReplaceMethod("B", "CNSet", function(object, value) {
	open(assayDataElement(object, "alleleA")) ## should do nothing if matrix
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

setMethod("flags", signature(object="CNSet"), function(object) flags(lM(object)))
setReplaceMethod("flags", signature=signature(object="CNSet", value="matrix"),
		 function(object, value){
			 linearParamElementReplace(object, "flags", value)
})

setAs("CNSetLM", "CNSet", function(from){
	if("batch" %in% varLabels(protocolData(from))){
		btch <- as.factor(protocolData(from)$batch)
	} else {
		stop("couldn't find batch in varLabels of protocolData.")
	}
	lm <- from@lM
	is.ffdf <- is(lm, "ffdf")
	if(is.ffdf){
		lm <- physical(lm)
	}
	lm.names <- c("tau2A", "tau2B", "sig2A", "sig2B", "nuA", "nuB", "phiA", "phiB", "phiPrimeA", "phiPrimeB", "corrAB", "corrAA", "corrBB")
	if(!all(lm.names %in% names(lm))){
		lm.names <- paste(lm.names, collapse=", ")
		stop(paste("names(object@lM) must have the following names:", lm.names))
	}
	tmp <- assayDataNew(tau2A=lm[["tau2A"]],
				tau2B=lm[["tau2B"]],
				sig2A=lm[["sig2A"]],
				sig2B=lm[["sig2B"]],
				nuA=lm[["nuA"]],
				nuB=lm[["nuB"]],
				phiA=lm[["phiA"]],
				phiB=lm[["phiB"]],
				phiPrimeA=lm[["phiPrimeA"]],
				phiPrimeB=lm[["phiPrimeB"]],
				corrAB=lm[["corrAB"]],
				corrAA=lm[["corrAA"]],
				corrBB=lm[["corrBB"]],
				flags=initializeBigMatrix("flags", nrow(from), length(unique(btch))))
	obj <- new("CNSet",
		   alleleA=assayData(from)[["alleleA"]],
		   alleleB=assayData(from)[["alleleB"]],
		   call=assayData(from)[["call"]],
		   callProbability=assayData(from)[["callProbability"]],
		   featureData=featureData(from),
		   phenoData=phenoData(from),
		   experimentData=experimentData(from),
		   protocolData=protocolData(from),
		   batch=btch,
		   lM=tmp)
	return(obj)
})
