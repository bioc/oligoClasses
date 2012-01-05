setMethod("isSnp", signature(object="gSet"),
	  function(object, ...) {
		  isSnp(featureData(object))
	  })

setMethod("isSnp", signature(object="character"),
	  function(object, pkgname, ...){
		  path <- system.file("extdata", package=pkgname)
		  load(file.path(path, "snpProbes.rda"))
		  snpProbes <- get("snpProbes")
		  res <- as.integer(object %in% snpProbes)
		  return(res)
	  })

setMethod("db", "gSet",
          function(object) {
		  requireAnnotation(annotation(object)) || stop(paste(annotation(object), "package not available"))
		  get(annotation(object))@getdb()
	  })

setMethod("chromosome", "gSet",
	  function(object, na.rm=FALSE, ...){
		  chromosome(featureData(object), na.rm)
	  })

setReplaceMethod("chromosome", signature(object="gSet", value="integer"),
		 function(object, value){
			 fData(object)$chromosome <-  value
			 object
		 })


setMethod("position", "gSet",
          function(object, na.rm=FALSE, ...){
		  position(featureData(object), na.rm)
          })

setMethod("checkOrder", signature(object="gSet"),
	  function(object, verbose=FALSE){
		  .checkOrder(object, verbose)
	  })


.checkOrder <- function(object, verbose=FALSE){
	d <- diff(order(chromosome(object), position(object)))
	if(any(d < 0)){
		if(verbose)
			warning("Object should be ordered by chromosome and physical position.\n",
				"Try \n",
				"> object <- order(object) \n")
		return(FALSE)
	}
	TRUE
}

chromosomePositionOrder <- function(object, ...){
	is.ordered <- checkOrder(object)
	if(!is.ordered){
		##if(verbose) message("Ordering ", class(object), " object by chromosome and physical position")
		index <- order(chromosome(object), position(object), ...)
		object <- object[index, ]
	}
	return(object)
}
