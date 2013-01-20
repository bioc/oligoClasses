setMethod("updateObject", signature(object="oligoSnpSet"),
          function(object, ..., verbose=FALSE) {
		  if (verbose) message("updateObject(object = 'oligoSnpSet')")
		  obj <- tryCatch(callNextMethod(object), error=function(e) NULL)
		  if(is.null(obj)){
			  obj <- new("oligoSnpSet",
				     assayData = updateObject(assayData(object),
				     ..., verbose=verbose),
				     phenoData = phenoData(object),
				     experimentData = updateObject(experimentData(object),
				     ..., verbose=verbose),
				     annotation = updateObject(annotation(object),
				     ..., verbose=verbose),
				     featureData=updateObject(featureData(object), ...,
				     genome=genomeBuild(object),
				     verbose=FALSE))
		  }
		  if (all(isCurrent(obj))) return(obj)
		  obj
          })

setReplaceMethod("baf", signature(object="oligoSnpSet"),
		     function(object, value){
			     if("baf" %in% assayDataElementNames(object)){
				     object <- assayDataElementReplace(object, "baf", value)
			     } else{
				     warning("assay data element 'baf' does not exist in oligoSnpSet object")
			     }
			     return(object)
		     })

setMethod("copyNumber", "oligoSnpSet", function(object) {
	cn <- assayDataElement(object, "copyNumber")
	return(cn)
})


setReplaceMethod("copyNumber", signature(object="oligoSnpSet", value="matrix"),
                 function(object, value){
			 assayDataElementReplace(object, "copyNumber", value)
		 })

setMethod("cnConfidence", "oligoSnpSet", function(object) assayData(object)[["cnConfidence"]])
setReplaceMethod("cnConfidence", signature(object="oligoSnpSet", value="matrix"),
                 function(object, value){
			 assayDataElementReplace(object, "cnConfidence", value)
                 })

setAs("oligoSnpSet", "data.frame",
      function(from, to){
	      cn <- copyNumber(from)/100
	      gt <- calls(from)
	      cn <- as.numeric(cn)
	      gt <- as.integer(gt)
	      baf.present <- "baf" %in% ls(assayData(from))
	      lrr.present <- "lrr" %in% ls(assayData(from))
	      if(baf.present){
		      bf <- as.numeric(assayDataElement(from, "baf"))/1000
	      }
	      if(lrr.present){
		      logRRatio <- as.numeric(assayDataElement(from, "lrr"))/100
	      }
	      x <- rep(position(from)/1e6, ncol(from))
	      ##x <- rep(position(object)[marker.index], 4)/1e6
	      is.snp <- rep(isSnp(from), ncol(from))
	      id <- rep(sampleNames(from), each=nrow(from))
	      if(!baf.present){
		      df <- data.frame(x=x, cn=cn, gt=gt, id=id,
				       is.snp=is.snp,
				       stringsAsFactors=FALSE)
	      } else {
		      df <- data.frame(x=x, cn=cn, gt=gt, baf=bf, id=id,
				       is.snp=is.snp,
				       stringsAsFactors=FALSE)
	      }
	      if(lrr.present){
		      df$lrr <- logRRatio
	      }
	      df$id <- factor(df$id, ordered=TRUE, levels=unique(df$id))
	      return(df)
      })

##setAs("oligoSnpSet", "SnpSet2", function(from, to){
##	new("SnpSet2",
##	    call=calls(from),
##	    callProbability=snpCallProbability(from),
##	    genome=genomeBuild(from),
##	    phenoData=phenoData(from),
##	    protocolData=protocolData(from),
##	    annotation=annotation(from),
##	    experimentData=experimentData(from),
##	    featureData=featureData(from))
##
##})

## ideally, oligoSnpSet would inherit from SnpSet2, but currently problems with GenomeAnnotatedDataFrame

##setMethod(snpCall, "oligoSnpSet", function(object, ...) {
##    assayDataElement(object, "call")
##})
##
##setMethod(snpCallProbability, "oligoSnpSet", function(object, ...) {
##    assayDataElement(object, "callProbability")
##})
##
##setReplaceMethod("snpCall", c("oligoSnpSet", "matrix"),
##                 function(object, ..., value)
##{
##    assayDataElementReplace(object, "call", value)
##})
##
##setReplaceMethod("snpCallProbability", c("oligoSnpSet", "matrix"),
##                 function(object, ..., value)
##{
##    assayDataElementReplace(object, "callProbability", value)
##})
##
####-----------------------
#### new methods for oligoSnpSet
####
##
##setMethod("calls", "oligoSnpSet", function(object) assayData(object)$call)
##setReplaceMethod("calls", signature(object="oligoSnpSet", value="matrix"),
##                 function(object, value)
##                 assayDataElementReplace(object, "call", value))
