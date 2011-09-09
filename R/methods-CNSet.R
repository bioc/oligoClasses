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
	x <- callNextMethod(x, i, j, ..., drop=FALSE)
	## ensure that assayData elements are matrices after subset operation
	isdf <- sapply(assayData(x), function(x) is(x, "data.frame"))
	if(any(isdf)){
		orig <- assayData(x)
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
	## one problem with the above -- the elements of assayData can be data.frame instead of matrix
##	phenoData(x) <- phenoData(x)[j, ...]
##	featureData(x) <- featureData(x)[i, ...]
##	protocolData(x) <- protocolData(x)[j, ...]
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




##setMethod("lM", "CNSet", function(object){
####	bs <- batchStatistics(object)
####	elem.names <- ls(bs)
####	param.names <- c("tau2A.AA", "tau2A.BB", "tau2B.AA", "tau2B.BB",
####			 "nuA", "nuB", "phiA", "phiB",
####			 "phiPrimeA", "phiPrimeB", "corrAB", "corrBB", "corrAA")
####	index <- match(param.names, elem.names)
####	return(bs[index])
##})
##setReplaceMethod("lM", signature=signature(object="CNSet", value="AssayData"),
##		 function(object, value){
##			 object@lM <- value
##			 object
##		 })

setMethod("batchStatistics", signature=signature(object="CNSet"), function(object) object@batchStatistics)
setReplaceMethod("batchStatistics", signature=signature(object="CNSet", value="AssayData"),
	 function(object, value){
		 object@batchStatistics <- value
		 object
	 })

##setMethod("Ns", signature=signature(object="CNSet"),
##	  function(object, batchname){
##		  if(missing(batchnames)) {
##			  stop("must specify batchname")
##		  } else {
##			  if(!batchname %in% batchNames(object))
##				  stop(paste("element.name must be one of ",  batchNames(object)))
##			  Ns <- matrix(NA, nrow(object), 3)
##			  colnames(NS) <- c("AA", "AB", "BB")
##			  j <- match(batchname, batchNames(object))
##			  Ns[, 1] <- assayDataElement(batchStatistics(object), "N.AA")[, j]
##			  Ns[, 2] <- assayDataElement(batchStatistics(object), "N.AB")[, j]
##			  Ns[, 3] <- assayDataElement(batchStatistics(object), "N.BB")[, j]
##			  return(Ns)
##		  }
##	  })
##setMethod("medians", signature=signature(object="CNSet"),
##	  function(object, batchname){
##		  if(missing(batchnames)) {
##			  stop("must specify batchname")
##		  } else {
##			  if(!batchname %in% batchNames(object))
##				  stop(paste("element.name must be one of ",  batchNames(object)))
##			  Ns <- matrix(NA, nrow(object), 3)
##			  colnames(NS) <- c("AA", "AB", "BB")
##			  j <- match(batchname, batchNames(object))
##			  Ns[, 1] <- assayDataElement(batchStatistics(object), "median.AA")[, j]
##			  Ns[, 2] <- assayDataElement(batchStatistics(object), "median.AB")[, j]
##			  Ns[, 3] <- assayDataElement(batchStatistics(object), "median.BB")[, j]
##			  return(Ns)
##		  }
##	  })
##setMethod("mads", signature=signature(object="CNSet"),
##	  function(object, batchname){
##		  if(missing(batchnames)) {
##			  stop("must specify batchname")
##		  } else {
##			  if(!batchname %in% batchNames(object))
##				  stop(paste("element.name must be one of ",  batchNames(object)))
##			  Ns <- matrix(NA, nrow(object), 3)
##			  colnames(NS) <- c("AA", "AB", "BB")
##			  j <- match(batchname, batchNames(object))
##			  Ns[, 1] <- assayDataElement(batchStatistics(object), "mad.AA")[, j]
##			  Ns[, 2] <- assayDataElement(batchStatistics(object), "mad.AB")[, j]
##			  Ns[, 3] <- assayDataElement(batchStatistics(object), "mad.BB")[, j]
##			  return(Ns)
##		  }
##	  })

setMethod("nu", c("CNSet", "character"), function(object, allele) nu(batchStatistics(object), allele))
setMethod("phi", c("CNSet", "character"), function(object, allele) phi(batchStatistics(object), allele))
##setMethod("phiPrime", c("CNSet", "character"), function(object, allele) phiPrime(batchStatistics(object), allele))
setMethod("sigma2", c("CNSet", "character"), function(object, allele) sigma2(batchStatistics(object), allele))
##setMethod("tau2", c("CNSet", "character"), function(object, allele) tau2(batchStatistics(object), allele))
##setMethod("corr", c("CNSet", "character"), function(object, allele) corr(batchStatistics(object), allele))
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
	##initialize container for Ns, but do not populate.
##	nr <- sum(fData(from)$isSnp, na.rm=TRUE)
##	Ns <- initializeNumberGenotypeFrom(from, unique(as.character(btch)), nr)
##	Ns <- initializeNumberGenotypeFrom(from, unique(as.character(btch)))
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

relocateObject <- function(object, to){
	stopifnot(isPackageLoaded("ff"))
	is.ff <- function(x) class(x)[1] == "ff_matrix"
	is.ffdf <- function(x) is(x, "ffdf")
	##AssayData
	storage.mode <- storageMode(assayData(object))
	orig <- assayData(object)
	physical <- get("physical")
	filename <- get("filename")
	pattern <- get("pattern")
	f2 <- function(X, dirname){
		##browser()
		X$filename <- file.path(dirname, basename(filename(X)))
		X$pattern <- file.path(dirname, basename(pattern(X)))
		physical(X)$filename <- file.path(dirname, basename(filename(X)))
		physical(X)$pattern <- file.path(dirname, basename(pattern(X)))
		X
	}
	assayData(object) <-
		switch(storage.mode,
		       environment=,
		       lockedEnvironment={
			       aData <- new.env(parent=emptyenv())
			       for (nm in ls(orig)) {
				       obj <- orig[[nm]]
				       if (is.ff(obj)) {
					       physical(obj)$pattern <- file.path(to, basename(pattern(obj)))
					       physical(obj)$filename = file.path(to, basename(filename(obj)))
				       }
				       else if (is.ffdf(obj)) {

					       ##important to not assign this to obj.  A list is returned (not a ffdf object, but the pointer is changed)
					       lapply(physical(obj), f2, dirname=to)
				       }
				       else message("Unable to change filename and pattern for assayData.")
				       aData[[nm]] <- obj
			       }
			       if ("lockedEnvironment" == storage.mode) Biobase:::assayDataEnvLock(aData)
			       aData
		       },
					#haven't tested yet when storage.mode returns list.
		       list = {
			       relocate.fxn <- function(obj,to){
				       obj <- orig[[nmm]]
				       filename(obj) <- file.path(to, basename(filename(obj)))
				       physical(obj)$pattern <- file.path(to, basename(pattern(obj)))
				       return(obj)
			       }
			       lapply(orig, relocate.fxn, to=to)
		       }
		       ) # ends switch
	storage.mode <- storage.mode(batchStatistics(object))
	orig <- batchStatistics(object)
	batchStatistics(object) <-
		switch(storage.mode,
		       environment =,
		       lockedEnvironment = {
			       aData <- new.env(parent=emptyenv())
			       for(nm in ls(orig)){
				       obj <- orig[[nm]]
				       if (is.ff(obj)){
					       physical(obj)$pattern = file.path(to, basename(pattern(obj)))
					       physical(obj)$filename = file.path(to, basename(filename(obj)))
				       }
				       else if (is.ffdf(obj)) {
					       f2 <- function(X, dirname){
						       X$filename <- file.path(dirname, basename(filename(X)))
						       X$pattern <- file.path(dirname, basename(pattern(X)))
						       physical(X)$filename <- file.path(dirname, basename(filename(X)))
						       physical(X)$pattern <- file.path(dirname, basename(pattern(X)))
						       X
					       }
					       ##important to not assign this to obj.  A list is returned (not a ffdf object, but the pointer is changed)
					       lapply(physical(obj), f2, dirname=to)
##					       fff <- function(X,dirname){
##						       physical(X)$pattern = file.path(dirname, basename(physical(X)$pattern))
##						       physical(X)$filename = file.path(dirname, basename(physical(X)$filename))
##						       return(X)
##					       }
##					       obj2 = lapply(obj, fff, dirname=to)
                                 }
                                 else message("Unable to change filename and pattern for batchStatistics.")
                                 aData[[nm]] <- obj
                             }
                             if ("lockedEnvironment" == storage.mode) Biobase:::assayDataEnvLock(aData)
                             aData
                         },
                         ## haven't tested yet when storage.mode returns list
                         list = {
                             relocate.fxn <- function(obj, to){
                                 obj <- orig[[nm]]
                                 filename(obj) <- file.path(to, basename(filename(obj)))
                                 physical(obj)$pattern <- file.path(to, basename(pattern(obj)))
                                 return(obj)
                             }
                             lapply(orig, relocate.fxn, to=to)
                         }
                         ) # ends switch
              ##Finally, remove any ff objects that might be in the phenodata
              tmp <- object$SNR
              if(is.ff(tmp)){
                  physical(tmp)$filename <- file.path(to, basename(filename(tmp)))
                  physical(tmp)$pattern <- file.path(to, basename(pattern(tmp)))
                  SNR <- tmp
                  tmp <- object$SKW
                  physical(tmp)$filename <- file.path(to, basename(filename(tmp)))
                  physical(tmp)$pattern <- file.path(to, basename(pattern(tmp)))
                  SKW <- tmp
                  pD <- new("AnnotatedDataFrame", data=data.frame(list(SNR=SNR[,],
                                                  SKW=SKW[,],
                                                  gender=object$gender)),
                            varMetadata=data.frame(labelDescription=c("SNR", "SKW", "gender"), row.names=c("SNR", "SKW", "gender")))
                  sampleNames(pD) <- sampleNames(object)
                  phenoData(object) <- pD
              }
              object
}




