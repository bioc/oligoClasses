setMethod("calls", "SnpSet", function(object) assayData(object)$call)
setReplaceMethod("calls", signature(object="SnpSet", value="matrix"),
                 function(object, value)
                 assayDataElementReplace(object, "call", value))

##setReplaceMethod("confs", signature(object="SnpSet", value="matrix"), function(object, value) assayDataElementReplace(object, "callProbability", value))

p2i <- function(p)
  as.integer(-1000*log(1-p))

i2p <- function(i)
  1-exp(-i/1000)

setMethod("confs", "SnpSet", function(object, transform=TRUE) {
	X <- assayData(object)$callProbability
	if (transform){
		X <- i2p(X)
	}
	return(X)
})

setReplaceMethod("confs", signature(object="SnpSet", value="matrix"),
		 function(object, value){
			 ##convert probability to integer
                         X <- matrix(p2i(value), nrow(X), ncol(X),
                                     dimnames=dimnames(value))
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
