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

##setMethod("callsConfidence", "SnpSet", function(object) confs(object))
##setReplaceMethod("callsConfidence", signature(object="SnpSet", value="matrix"),
##                 function(object, value) confs(object) <- value)

setMethod("combine", signature=signature(x="SnpSet", y="SnpSet"),
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

## need this to work for a RangedData object with multiple ranges
setMethod("featuresInRange", signature(object="SnpSet", range="RangedDataCNV"),
	  function(object, range, FRAME=0, FRAME.LEFT, FRAME.RIGHT, ...){
		  start <- start(range)
		  end <- end(range)
		  CHR <- chromosome(range)
		  ##featuresInXlim(object, start=start(range), end=end(range), CHR=range$chrom, ...)
		  if(missing(FRAME.LEFT)) FRAME.LEFT <- FRAME
		  if(missing(FRAME.RIGHT)) FRAME.RIGHT <- FRAME
		  data(chromosomeAnnotation, package="SNPchip")
		  chr.end <- chromosomeAnnotation[CHR, "chromosomeSize"]
		  start <- max(start-FRAME.LEFT, 0)
		  end <- min(end+FRAME.RIGHT, chr.end)
		  which(position(object) >= start & position(object) <= end & chromosome(object) == CHR)
	  })

##setMethod("order", "SnpSet",
##	  function(..., na.last=TRUE, decreasing=FALSE){
##		  chromosomePositionOrder(...)
##	  })


setMethod("checkOrder", signature(object="SnpSet"),
	  function(object, verbose=FALSE){
		  .checkOrder(object, verbose)
	  })
