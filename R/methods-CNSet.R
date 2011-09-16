setMethod("show", "CNSet", function(object){
	is.ff <- is(calls(object), "ff_matrix") | is(calls(object), "ffdf")
	if(is.ff){
		##to avoid warnings
		if("SKW" %in% varLabels(object)) {
			if(is(object$SKW, "ff"))
				open(object$SKW)
		}
		if("SNR" %in% varLabels(object)){
			if(is(object$SNR, "ff"))
				open(object$SNR)
		}
		if("gender" %in% varLabels(object)){
			if(is(object$gender, "ff"))
				open(object$gender)
		}
	}
	ad.class <- class(A(object))[1]
	cat("CNSet (assayData/batchStatistics elements: ", ad.class, ")\n", sep="")
	callNextMethod(object)
	bns <- batchNames(object)
	bns <- bns[-length(bns)]
	freq <- as.integer(table(batch(object)))
	cat("batch:   ", paste(bns, freq, sep=", "), "\n")
	adim <- list(nrow(object), length(batchNames(object)))
	cat("batchStatistics: ", length(ls(batchStatistics(object))), " elements, ", nrow(object), " features, ", length(unique(batch(object))), " batches\n")
})

setMethod("[", "CNSet", function(x, i, j, ..., drop=FALSE){
	##isff <- is(A(x), "ff") | is(A(x), "ffdf")
	##if(isff) open(x)
	open(x)
	x <- callNextMethod(x, i, j, ..., drop=FALSE)
##	if(isff) close(x)
##	x <- xx; rm(x)
	## ensure that assayData elements are matrices after subset operation
	isdf <- is(A(x), "data.frame")
	if(isdf){
		orig <- assayData(x)
		storage.mode <- Biobase:::assayDataStorageMode(orig)
		assayData(x) <-
			switch(storage.mode,
			       environment =,
			       lockedEnvironment = {
				       aData <- new.env(parent=emptyenv())
				       for(nm in ls(orig)) aData[[nm]] <- as.matrix(orig[[nm]])##[i, j, ..., drop = drop]
				       if ("lockedEnvironment" == storage.mode) Biobase:::assayDataEnvLock(aData)
				       aData
			       },
			       list = {
				       lapply(orig, as.matrix)
			       })
	}
	if(missing(j)) j <- 1:ncol(x)
	if(missing(i)) i <- 1:nrow(x)
	x@batch <- batch(x)[j]
	nms <- sampleNames(batchStatistics(x))
	## need to subset columns of LinearModelParameter
	## Adapted from the '[' method for eSet in Biobase
	## redefine 'j'
	j <- which(nms %in% unique(as.character(batch(x))))
	storage.mode <- storageMode(batchStatistics(x))
	## i (if defined) is already subset by callNextMethod
	orig <- batchStatistics(x)
	batchStatistics(x) <-
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
	return(x)
})

setMethod("batch", "CNSet", function(object) object@batch)

setReplaceMethod("batch", signature=signature(object="CNSet"),
	 function(object, value){
		 object@batch <- as.character(value)
		 object
})


setMethod("batchNames", "CNSet", function(object)  batchNames(batchStatistics(object)))

setReplaceMethod("batchNames", "CNSet", function(object, value) {
	batchNames(batchStatistics(object)) <- value
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
	obj <- assayDataElementReplace(object, "alleleA", value)
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
	names <- ls(batchStatistics(object))
	L <- length(names)
	for(i in 1:L) {
		tmp <- eval(substitute(assayData(object)[[NAME]], list(NAME=names[i])))
		if(!is.null(tmp)) close(tmp)
	}
	##lapply(physical(batchStatistics(con)), open)
	return()
})

setMethod("open", "CNSet", function(con, ...){
	object <- con
	if(!isFF(object)) return()
	names <- ls(assayData(object))
	L <- length(names)
	for(i in 1:L) open(eval(substitute(assayData(object)[[NAME]], list(NAME=names[i]))))
	names <- assayDataElementNames(batchStatistics(object))
	L <- length(names)
	for(i in 1:L) open(eval(substitute(batchStatistics(object)[[NAME]], list(NAME=names[i]))))
	if("SKW" %in% varLabels(object)){
		if(is(object$SKW, "ff")) open(object$SKW)
	}
	if("SNR" %in% varLabels(object)){
		if(is(object$SNR, "ff")) open(object$SNR)
	}
	return(TRUE)
})
## does nothing if not an ff object
setMethod("open", "numeric", function(con, ...) return(con))
setMethod("open", "matrix", function(con, ...) return(con))
setMethod("open", "array", function(con, ...) return(con))
setMethod("close", "numeric", function(con, ...) return(con))
setMethod("close", "matrix", function(con, ...) return(con))
setMethod("close", "array", function(con, ...) return(con))



setMethod("nu", c("CNSet", "character"), function(object, allele) nu(batchStatistics(object), allele))
setMethod("phi", c("CNSet", "character"), function(object, allele) phi(batchStatistics(object), allele))
setMethod("sigma2", c("CNSet", "character"), function(object, allele) sigma2(batchStatistics(object), allele))
setMethod("flags", signature(object="CNSet"), function(object) flags(batchStatistics(object)))

setAs("CNSetLM", "CNSet", function(from){
	if("batch" %in% varLabels(protocolData(from))){
		btch <- as.character(protocolData(from)$batch)
	} else {
		stop("couldn't find batch in varLabels of protocolData.")
	}
	lm <- from@lM
	is.ffdf <- is(lm, "ffdf")
	if(is.ffdf){
		##stopifnot(isPackageLoaded("ff"))
		lm <- physical(lm)
	}
	lm.names <- c("tau2A", "tau2B", "sig2A", "sig2B", "nuA", "nuB", "phiA", "phiB", "phiPrimeA", "phiPrimeB", "corrAB", "corrAA", "corrBB")
	if(!all(lm.names %in% names(lm))){
		lm.names <- paste(lm.names, collapse=", ")
		stop(paste("names(object@lM) must have the following names:", lm.names))
	}
	nr <- nrow(from)
	nc <- length(unique(btch))
	## mainly to avoid initializing new ff objects
	tau2A.AA <- lm[["sig2A"]]
	tau2A.BB <- lm[["tau2A"]]
	tau2B.AA <- lm[["tau2B"]]
	tau2B.BB <- lm[["sig2B"]]
	tmp <- assayDataNew(N.AA=initializeBigMatrix("N.AA", nr, nc),
			    N.AB=initializeBigMatrix("N.AB", nr, nc),
			    N.BB=initializeBigMatrix("N.BB", nr, nc),
			    medianA.AA=initializeBigMatrix("median.AA", nr, nc),
			    medianA.AB=initializeBigMatrix("median.AB", nr, nc),
			    medianA.BB=initializeBigMatrix("median.BB", nr, nc),
			    medianB.AA=initializeBigMatrix("median.AA", nr, nc),
			    medianB.AB=initializeBigMatrix("median.AB", nr, nc),
			    medianB.BB=initializeBigMatrix("median.BB", nr, nc),
			    madA.AA=initializeBigMatrix("mad.AA", nr, nc),
			    madA.AB=initializeBigMatrix("mad.AB", nr, nc),
			    madA.BB=initializeBigMatrix("mad.BB", nr, nc),
			    madB.AA=initializeBigMatrix("mad.AA", nr, nc),
			    madB.AB=initializeBigMatrix("mad.AB", nr, nc),
			    madB.BB=initializeBigMatrix("mad.BB", nr, nc),
			    tau2A.AA=tau2A.AA,
			    tau2A.BB=tau2A.BB,
			    tau2B.AA=tau2B.AA,
			    tau2B.BB=tau2B.BB,
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
		   batchStatistics=tmp)
	return(obj)
})

setMethod("batchStatistics", signature=signature(object="CNSet"), function(object) object@batchStatistics)
setReplaceMethod("batchStatistics", signature=signature(object="CNSet", value="AssayData"),
	 function(object, value){
		 object@batchStatistics <- value
		 object
	 })





