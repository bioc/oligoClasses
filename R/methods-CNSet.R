##setMethod("CA", "CNSet", function(object) assayData(object)[["CA"]]/100)
##setMethod("CB", "CNSet", function(object) assayData(object)[["CB"]]/100)
setMethod("CB", "CNSet", function(object){
	assayDataElement(object, "CB")
})
setMethod("CA", "CNSet", function(object) {
	assayDataElement(object, "CA")
})
setReplaceMethod("CB", "CNSet", function(object, value) assayDataElementReplace(object, "CB", value))
setReplaceMethod("CA", "CNSet", function(object, value) assayDataElementReplace(object, "CA", value))
##setReplaceMethod("CA", signature(object="CNSet", value="matrix"),
##		 function(object, value){
##			 value <- matrix(as.integer(value*100), nrow(value), ncol(value), dimnames=dimnames(value))
##			 assayDataElementReplace(object, "CA", value)		
##		 })
##setReplaceMethod("CB", signature(object="CNSet", value="matrix"),
##		 function(object, value){
##			 value <- matrix(as.integer(value*100), nrow(value), ncol(value), dimnames=dimnames(value))			 
##			 assayDataElementReplace(object, "CB", value)
##		 })

setMethod("totalCopyNumber", "CNSet", function(object, i, j){
	if(missing(i) & missing(j)){
		if(inherits(CA(object), "ff") | inherits(CA(object), "ffdf")) stop("Must specify i and/or j for ff objects")
	}
	if(missing(i) & !missing(j)){
		snp.index <- which(isSnp(object))	
		cn.total <- as.matrix(CA(cnSet)[, j])
		cb <- as.matrix(CB(cnSet)[snp.index, j]	)
		cn.total[snp.index, ] <- cn.total[snp.index, ] + cb		
	}
	if(!missing(i) & missing(j)){
		snp.index <- intersect(which(isSnp(object)), i)
		cn.total <- as.matrix(CA(cnSet)[i, ])
		cb <- as.matrix(CB(cnSet)[snp.index, ])	
		cn.total[snp.index, ] <- cn.total[snp.index, ] + cb				
	}
	if(!missing(i) & !missing(j)){
		snp.index <- intersect(which(isSnp(object)), i)		
		cn.total <- as.matrix(CA(cnSet)[i, j])	
		cb <- as.matrix(CB(cnSet)[snp.index, j])
		cn.total[snp.index, ] <- cn.total[snp.index, ] + cb
	}
	cn.total <- cn.total/100
	dimnames(cn.total) <- NULL
	return(cn.total)
})


