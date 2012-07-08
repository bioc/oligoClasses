##
## Directly from Biobase
##
setMethod("initialize", "SnpSet2",
          function(.Object,
                   assayData = assayDataNew(call = call,
                                            callProbability = callProbability, ...),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData,## = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   protocolData = phenoData[,integer(0)],
                   call = new("matrix"),
                   callProbability = matrix(numeric(),
                                            nrow=nrow(call), ncol=ncol(call),
                                            dimnames=dimnames(call)),
		   genome=c("hg19", "hg18"),
                   ...) {
		  genome <- match.arg(genome)
		  if(missing(featureData))
			  featureData <- GenomeAnnotatedDataFrameFrom(assayData, annotation, genome=genome)
		  callNextMethod(.Object,
				 assayData = assayData,
				 phenoData = phenoData,
				 featureData = featureData,
				 experimentData = experimentData,
				 annotation = annotation,
				 protocolData = protocolData, ...)
          })

setMethod(snpCall, "SnpSet2", function(object, ...) {
    assayDataElement(object, "call")
})

setMethod(snpCallProbability, "SnpSet2", function(object, ...) {
    assayDataElement(object, "callProbability")
})



setReplaceMethod("snpCall", c("SnpSet2", "matrix"),
                 function(object, ..., value){
			 assayDataElementReplace(object, "call", value)
		 })

setReplaceMethod("snpCallProbability", c("SnpSet2", "matrix"),
                 function(object, ..., value){
			 assayDataElementReplace(object, "callProbability", value)
		 })

##-----------------------
## new methods for SnpSet2
##
setMethod("calls", "SnpSet2", function(object) assayData(object)$call)
setReplaceMethod("calls", signature(object="SnpSet2", value="matrix"),
                 function(object, value)
                 assayDataElementReplace(object, "call", value))

setMethod("calls", "SnpSet", function(object) assayData(object)$call)
setReplaceMethod("calls", signature(object="SnpSet", value="matrix"),
                 function(object, value)
                 assayDataElementReplace(object, "call", value))


p2i <- function(p)
  as.integer(-1000*log(1-p))

i2p <- function(i)
  1-exp(-i/1000)

warningMsg <- function(X){
	.class=class(X)
	warning("callProbability slot is of class ", .class, ".\n")
	cat("\nTo obtain the confidence scores, the data needs to be extracted from disk and represented as a matrix. The '[' method does both.  For example,\n", fill=TRUE)
	message("> x <- confs(object)[,] ## 'x' is a matrix\n")
	cat("* Note however that 'x' may be very large and swamp  the available RAM. A better approach would be to specify which rows (i) and columns (j) are read only those rows and columns from disk.\n", fill=TRUE)
	message("> x < confs(object)[i, j] \n")
	message("Finally, 'x' still needs to be translated to a probability.  This can be done by", fill=TRUE)
	message("> p <- i2p(x)")
}

setMethod("confs", "SnpSet2", function(object, transform=TRUE) {
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

setReplaceMethod("confs", signature(object="SnpSet2", value="matrix"),
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

setMethod("confs", "SnpSet", function(object, transform=TRUE) {
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

setReplaceMethod("confs", signature(object="SnpSet", value="matrix"),
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



setMethod("combine", signature=signature(x="SnpSet2", y="SnpSet2"),
          function(x, y, ...){
		  ##Check that both x and y are valid objects
		  if(!validObject(x)) stop("x is not a valid object")
		  if(!validObject(y)) stop("y is not a valid object")
		  annot <- paste(sort(c(annotation(x), annotation(y))), collapse=",")
		  annotation(x) <- annotation(y) <- annot

		  if(class(x) != class(y)){
			  stop("objects must have the same class")
		  }
		  if(storageMode(assayData(x)) != storageMode(assayData(y))){
			  stop("objects must have same storage mode for assayData")
		  }

		  fd <- combine(featureData(x), featureData(y))
		  pd <- combine(phenoData(x), phenoData(y))
		  ad.x <- as.list(assayData(x))
		  ad.y <- as.list(assayData(y))
		  ad.xy <- mapply(rbind, ad.x, ad.y, SIMPLIFY=FALSE)
		  id.x <- match(rownames(ad.xy[[1]]), featureNames(fd))
		  ee <- combine(experimentData(x), experimentData(y))
		  assayData(x) <- ad.xy
		  storageMode(assayData(x)) <- storageMode(assayData(y))
		  experimentData(x) <- ee
		  featureData(x) <- fd
		  phenoData(x) <- pd
		  x
          })


setMethod("featuresInRange", signature(object="SnpSet2", range="RangedDataCNV"),
	  function(object, range, FRAME=0, FRAME.LEFT, FRAME.RIGHT, ...){
		  .Defunct("featuresInRange has been deprecated. Use findOverlaps.")
	  })

