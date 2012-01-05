setMethod("chromosome", signature(object="AnnotatedDataFrame"),
	  function(object, na.rm=FALSE, ...){
		  chrom <- object$chromosome
		  if(!na.rm) return(chrom)
		  chrom[!is.na(chrom)]
	  })

setMethod("position", signature(object="AnnotatedDataFrame"),
	  function(object, na.rm=FALSE, ...) {
		  pos <- object$position
		  if(!na.rm){
			  return(pos)
		  }
		  pos[!is.na(pos)]
	})
