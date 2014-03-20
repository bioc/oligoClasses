setMethod("initialize", signature(.Object="BeadStudioSetList"),
	  function(.Object,
		   assayDataList=AssayDataList(baf=baf, lrr=lrr),
		   lrr=list(),
		   baf=lapply(lrr, function(x) matrix(nrow=nrow(x), ncol=ncol(x))),
		   featureDataList=GenomeAnnotatedDataFrameFrom(assayDataList, annotation, genome),
		   chromosome=vector("list", length(lrr)),
		   phenoData,
		   annotation=character(),
		   genome=character(),
		   ...){
		  if(missing(phenoData)){
			  if(length(lrr) > 0){
				  phenoData <- annotatedDataFrameFrom(lrr[[1]], byrow=FALSE)
			  } else {
				  phenoData <- new("AnnotatedDataFrame")
			  }
		  }
		  callNextMethod(.Object,
				 assayDataList=assayDataList,
				 featureDataList=featureDataList,
				 phenoData=phenoData,
				 chromosome=chromosome,
				 annotation=annotation,
				 genome=genome,
				 ...)
	  })

setMethod("updateObject", signature(object="BeadStudioSetList"),
          function(object, ..., verbose=FALSE) {
		  if (verbose) message("updateObject(object = 'BeadStudioSetList')")
		  obj <- tryCatch(callNextMethod(object), error=function(e) NULL)
		  if(is.null(obj)){
			  obj <- new("BeadStudioSetList",
				     assayDataList = assayDataList(object),
				     phenoData = phenoData(object),
				     annotation = updateObject(annotation(object),
				     ..., verbose=verbose),
				     featureDataList=featureDataList(object),
				     chromosome=chromosome(object),
				     genome=genomeBuild(object),
				     ...)
		  }
		  obj
          })
setMethod("[[", signature(x="BeadStudioSetList"),
	  function(x, i, j, ..., exact=TRUE){
		  if(missing(i)) return(x)
		  ad <- assayDataList(x)
		  fdlist <- featureData(x)
		  adnew <- switch(storage.mode(ad),
			  lockedEnvironment =,
				  environment = new.env(parent=emptyenv()),
				  list = list())
		  nms <- ls(ad)
		  if(length(i) == 1){
			  for (nm in ls(ad)){
				  elt <- ad[[nm]][[i]]
				  dimnames(elt) <- lapply(dimnames(elt), unname)
				  adnew[[nm]] <- elt
			  }
		  }
		  x <- new("BeadStudioSet",
			   assayData=adnew,
			   phenoData=phenoData(x),
			   featureData=fdlist[[i]],
			   genome=genomeBuild(x),
			   annotation=annotation(x))
	  })

setMethod("[[", signature(x="BafLrrSetList"),
	  function(x, i, j, ..., exact=TRUE){
		 x <- callNextMethod()
		 new("BafLrrSet",
		     assayData=assayData(x),
		     phenoData=phenoData(x),
		     featureData=featureData(x),
		     genome=genomeBuild(x),
		     annotation=annotation(x))
	  })

setMethod("[", signature(x="gSetList"),
	  function(x, i, j, ..., drop=TRUE){
		  if(missing(i) && missing(j)) return(x)
		  ad <- assayDataList(x)
		  if(!missing(i)){
			  fdlist <- featureData(x)[i]
		  }
		  adnew <- switch(storage.mode(ad),
			  lockedEnvironment =,
				  environment = new.env(parent=emptyenv()),
				  list = list())
		  nms <- ls(ad)
		  if(!missing(i)){
			  for (nm in ls(ad)){
				  elt <- ad[[nm]][i]
				  adnew[[nm]] <- elt
			  }
			  ad <- adnew
		  }
		  if(missing(j)){
			  x@featureDataList <- fdlist
			  x@chromosome <- x@chromosome[i]
		  } else {
			  for (nm in ls(ad)){
				  elt <- lapply(ad[[nm]], function(y, j) y[, j, drop=FALSE], j=j)
				  adnew[[nm]] <- elt
			  }
			  phenoData(x) <- phenoData(x)[j, ]
			  x@protocolData <- x@protocolData[j, ]
		  }
		  x@assayDataList <- adnew
		  return(x)
	  })

setReplaceMethod("[[", signature(x="BafLrrSetList", value="BafLrrSet"),
		 function(x, i, j, ..., value){
			 fdl <- x@featureDataList
			 fdl[[i]] <- featureData(value)
			 adl <- x@assayDataList
			 r <- adl[["lrr"]]
			 r[[i]] <- lrr(value)
			 b <- adl[["baf"]]
			 b[[i]] <- baf(value)
			 adl <- AssayDataList(lrr=r, baf=b)
			 new("BafLrrSetList",
			     assayDataList=adl,
			     featureDataList=fdl,
			     phenoData=phenoData(x),
			     chromosome=chromosome(x),
			     annotation=annotation(x),
			     genome=genomeBuild(x))
		 })


setReplaceMethod("assayData", signature=signature(object="BeadStudioSetList",
			      value="AssayData"),
                 function(object, value) {
			 object@assayDataList <- value
			 object
                 })


setMethod("baf", signature(object="oligoSetList"),
	  function(object) assayData(object)[["baf"]])
setMethod("baf", signature(object="BeadStudioSetList"),
	  function(object){
		  ##lapply(object, baf)
		  assayDataList(object)[["baf"]]
	  })

setMethod(baf, signature(object="BafLrrSetList"),
	  function(object){
            assayDataList(object)[["baf"]]
	  })

setMethod("calls", signature(object="oligoSetList"),
	  function(object) assayData(object)[["call"]])
setMethod("copyNumber", signature(object="oligoSetList"),
	  function(object) assayData(object)[["copyNumber"]])

setMethod("lrr", signature(object="BeadStudioSetList"),
	  function(object){
		  ##lapply(object, lrr)
		  assayDataList(object)[["lrr"]]
	  })

setMethod("lrr", signature(object="BafLrrSetList"),
	  function(object){
            ##lapply(object, lrr)
            assayDataList(object)[["lrr"]]
	  })


setReplaceMethod("lrr", signature(object="BafLrrSetList", value="matrix"),
	  function(object, value){
		  ## value can often be fewer columns than object
		  if(is.null(rownames(value))) stop("row.names is NULL")
		  if(is.null(colnames(value))) stop("col.names is NULL")
		  sample.index <- match(colnames(value), sampleNames(object))
		  for(j in seq_along(object)){
			  bset <- object[[j]]
			  k <- match(featureNames(bset), rownames(value))
			  lrr(bset)[, sample.index] <- value[k, , drop=FALSE]
			  object[[j]] <- bset
		  }
		  return(object)
	  })

##setMethod("ncol", signature(x="BeadStudioSetList"),
##	  function(x) ncol(assayDataList(x)[["lrr"]][[1]]))

setMethod("snpCallProbability", signature(object="oligoSetList"),
	  function(object) assayData(object)[["callProbability"]])


setMethod(clone, "BafLrrSetList", function(object, id, prefix, ...){
	duplicateBLList(object, ids=id, prefix=prefix, ...)
})

duplicateBLList <- function(object, ids, prefix="waveAdj", empty=FALSE){
	##brList.copy <- object
	## duplicate the lrr ff objects.  Then do wave correction on the
	## duplicated files.
	if(missing(ids)) ids <- sampleNames(object)
	ids <- as.character(ids)
	r <- lrr(object)
	b <- baf(object)
	rcopy.list <- list()
	bcopy.list <- list()
	for(i in seq_along(r)){
		x <- r[[i]]
		y <- b[[i]]
		rcopy <- initializeBigMatrix(paste(prefix, "lrr", sep="-"), nrow(x), length(ids), vmode="integer")
		bcopy <- initializeBigMatrix(paste(prefix, "baf", sep="-"), nrow(x), length(ids), vmode="integer")
		dimnames(rcopy) <- list(rownames(x),
					ids)
		dimnames(bcopy) <- dimnames(rcopy)
		J <- match(ids, colnames(x))
		if(!empty){
			for(j in seq_along(J)){
				k <- J[j]
				rcopy[, j] <- x[, k]
				bcopy[, j] <- y[, k]
			}
		}
		rcopy.list[[i]] <- rcopy
		bcopy.list[[i]] <- bcopy
	}
	adl <- AssayDataList(baf=bcopy.list, lrr=rcopy.list)
	pd <- phenoData(object)[match(ids, sampleNames(object)), ]
	new("BafLrrSetList",
	    assayDataList=adl,
	    featureDataList=featureData(object),
	    phenoData=pd,
	    chromosome=chromosome(object),
	    annotation=annotation(object),
	    genome=genomeBuild(object))
}
