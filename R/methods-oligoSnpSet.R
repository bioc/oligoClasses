##setMethod("copyNumber", "oligoSnpSet", function(object) {
##	cn <- assayDataElement(object, "copyNumber")
####	if(is(cn, "numeric")) {
####		return(cn)
####	}
####	if(is(cn, "integer")){
####		cn <- cn/100
####		return(cn)
####	}
##	return(cn)
##})

setMethod("copyNumber", "oligoSnpSet", function(object) {
	cn <- assayDataElement(object, "copyNumber")
	return(cn)
})


setReplaceMethod("copyNumber", signature(object="oligoSnpSet", value="matrix"),
                 function(object, value){
			 ##value <- matrix(as.integer(value*100), nrow(value), ncol(value), dimnames=dimnames(value))
			 assayDataElementReplace(object, "copyNumber", value)
		 })

setMethod("cnConfidence", "oligoSnpSet", function(object) assayData(object)[["cnConfidence"]])
setReplaceMethod("cnConfidence", signature(object="oligoSnpSet", value="matrix"),
                 function(object, value){
##			 conf <- value
##			 dns <- dimnames(conf)
##			 X <- matrix(as.integer(conf*100), nrow(conf), ncol(conf))
##			 dimnames(X) <- dns
			 assayDataElementReplace(object, "cnConfidence", value)
                 })

setAs("oligoSnpSet", "data.frame",
      function(from, to){
	      cn <- copyNumber(from)
	      gt <- calls(from)
	      cn <- as.numeric(cn)
	      gt <- as.integer(gt)
	      baf.present <- "baf" %in% ls(assayData(from))
	      if(baf.present){
		      bf <- as.numeric(assayDataElement(from, "baf"))
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
	      return(df)
      })
