setMethod("chromosome", signature(object="SummarizedExperiment"),
	  function(object,...) as.character(seqnames(object)))
setMethod("isSnp", signature(object="SummarizedExperiment"),
	  function(object,...) values(rowData(object))$isSnp)
setMethod("lrr", signature(object="SummarizedExperiment"),
	  function(object){
		  assays(object)[[1]]
	  })
setMethod("baf", signature(object="SummarizedExperiment"),
	  function(object){
		  assays(object)[[2]]
	  })
