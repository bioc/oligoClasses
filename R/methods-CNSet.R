setMethod("show", "CNSet", function(object){
	is.ff <- is(calls(object), "ff_matrix") | is(calls(object), "ffdf")
	if(is.ff){
		if(!isPackageLoaded("ff")) warning("ff objects detected, but ff package is not loaded")
		##to avoid warnings
		if("SKW" %in% varLabels(object)) {
			if(is(object$SKW, "ff"))
				open.ff(object$SKW)
		}
		if("SNR" %in% varLabels(object)){
			if(is(object$SNR, "ff"))
				open.ff(object$SNR)
		}
		if("gender" %in% varLabels(object)){
			if(is(object$gender, "ff"))
				open.ff(object$gender)
		}
	}
	ad.class <- class(A(object))[1]
	cat("CNSet (assayData/batchStatistics elements: ", ad.class, ")\n", sep="")
	callNextMethod(object)
	bns <- head(batchNames(object))
	##bns <- bns[-length(bns)]
	freq <- as.integer(table(batch(object)))
	index <- names(table(batch(object))) %in% bns
	freq <- freq[index]
	cat("batch:   ", paste(bns, ":", freq, sep="", collapse=", "), "\n")
	adim <- list(nrow(object), length(batchNames(object)))
	if(adim[[1]] > 0){
		cat("batchStatistics: ", length(ls(batchStatistics(object))), " elements, ", nrow(object), " features, ", length(unique(batch(object))), " batches\n")
	}
})

setMethod("updateObject", signature(object="CNSet"),
          function(object, ..., verbose=FALSE) {
		  if (verbose) message("updateObject(object = 'CNSet')")
		  obj <- tryCatch(callNextMethod(batch=batch(object)), error=function(e) NULL)
		  if(is.null(obj)){
			  ## must supply batch for batchStatistics to be added
			  if(is(calls(object), "ffdf") | is(calls(object), "ff_matrix"))
				  stopifnot(isPackageLoaded("ff"))
			  if(.hasSlot(object, "mixtureParams")){
				  obj <- new("CNSet",
					     assayData = updateObject(assayData(object),
					     ..., verbose=verbose),
					     phenoData = phenoData(object),
					     experimentData = experimentData(object),
					     annotation = updateObject(annotation(object),
					     ..., verbose=verbose),
					     featureData=updateObject(featureData(object), ..., verbose=verbose),
					     batch=as.character(batch(object)),
					     batchStatistics=batchStatistics(object),
					     mixtureParams=object@mixtureParams)
			  } else {
				  obj <- new("CNSet",
					     assayData = updateObject(assayData(object),
					     ..., verbose=verbose),
					     phenoData = phenoData(object),
					     experimentData = experimentData(object),
					     annotation = updateObject(annotation(object),
					     ..., verbose=verbose),
					     featureData=updateObject(featureData(object), ..., verbose=verbose),
					     batch=as.character(batch(object)),
					     batchStatistics=batchStatistics(object),
					     mixtureParams=matrix(NA, 4, ncol(object)))
			  }
			  if (isCurrent(obj)["CNSet"]) return(obj)
			  return(obj)
		  }
          })

setMethod("[", "CNSet", function(x, i, j, ..., drop=FALSE){
  openff(x)
  x <- callNextMethod(x, i, j, ..., drop=FALSE)
  isdf <- is(A(x), "data.frame")
  if(isdf){
    orig <- assayData(x)
    ##storage.mode <- Biobase:::assayDataStorageMode(orig)
    storage.mode <- storageMode(orig)
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

setMethod("closeff",
	  signature(object="CNSet"), function(object){
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
})

setMethod("close", "CNSet", function(con, ...){
	object <- con
	closeff(object)
})

setMethod("openff", signature(object="CNSet"),
	  function(object){
		  if(isFF(object)){
			  names <- ls(assayData(object))
			  L <- length(names)
			  for(i in 1:L) open(eval(substitute(assayData(object)[[NAME]], list(NAME=names[i]))))
			  names <- assayDataElementNames(batchStatistics(object))
			  L <- length(names)
			  for(i in 1:L) open(eval(substitute(batchStatistics(object)[[NAME]], list(NAME=names[i]))))
		  }
		  if("SKW" %in% varLabels(object)){
			  if(is(object$SKW, "ff")) open(object$SKW)
		  }
		  if("SNR" %in% varLabels(object)){
			  if(is(object$SNR, "ff")) open(object$SNR)
		  }
		  if("SNR" %in% varLabels(object)){
			  if(is(object$gender, "ff")) open(object$gender)
		  }
	  })

setMethod("open", "CNSet", function(con, ...){
	object <- con
	openff(object)
	return(TRUE)
})

setMethod("nu", c("CNSet", "character"), function(object, allele) nu(batchStatistics(object), allele))
setMethod("phi", c("CNSet", "character"), function(object, allele) phi(batchStatistics(object), allele))
setMethod("sigma2", c("CNSet", "character"), function(object, allele) sigma2(batchStatistics(object), allele))
setMethod("flags", signature(object="CNSet"), function(object) flags(batchStatistics(object)))

setMethod("batchStatistics", signature=signature(object="CNSet"), function(object) object@batchStatistics)
setReplaceMethod("batchStatistics", signature=signature(object="CNSet", value="AssayData"),
	 function(object, value){
		 object@batchStatistics <- value
		 object
	 })


setMethod(snpCall, "CNSet", function(object, ...) {
	assayDataElement(object, "call")
})

setMethod(snpCallProbability, "CNSet", function(object, ...) {
	assayDataElement(object, "callProbability")
})

setReplaceMethod("snpCall", c("CNSet", "matrix"),
                 function(object, ..., value)
	 {
		 assayDataElementReplace(object, "call", value)
	 })

setReplaceMethod("snpCallProbability", c("CNSet", "matrix"),
                 function(object, ..., value)
{
	assayDataElementReplace(object, "callProbability", value)
})

setMethod("calls", "CNSet", function(object) assayData(object)$call)
setReplaceMethod("calls", signature(object="CNSet", value="matrix"),
                 function(object, value)
                 assayDataElementReplace(object, "call", value))

setMethod("confs", "CNSet", function(object, transform=TRUE) {
	X <- snpCallProbability(object)
	if(is(X, "ff_matrix") | is(X, "ffdf")){
		warningMsg(X)
		return(X)
	}
	if (transform){
		X <- i2p(X)
	}
	return(X)
})

setReplaceMethod("confs", signature(object="CNSet", value="matrix"),
		 function(object, value){
			 ##convert probability to integer
			 if(max(value) > 1){
				 X <- matrix(p2i(value), nrow(X), ncol(X),
					     dimnames=dimnames(value))
			 } else {
				 X <- value
			 }
			 assayDataElementReplace(object, "callProbability", X)
		 })

##setMethod("genomeBuild", "CNSet", function(object) object@genome)
##setReplaceMethod("genomeBuild", signature(object="CNSet", value="character"),
##		 function(object, value){
##			 object@genome <- value
##			 return(object)
##		 })
