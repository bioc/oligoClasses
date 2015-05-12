setMethod("chromosome", signature(object="RangedSummarizedExperiment"),
	  function(object,...) as.character(seqnames(object)))
setMethod("isSnp", signature(object="RangedSummarizedExperiment"),
	  function(object,...) values(rowRanges(object))$isSnp)
setMethod("lrr", signature(object="RangedSummarizedExperiment"),
	  function(object){
		  assays(object)[[1]]
	  })
setMethod("baf", signature(object="RangedSummarizedExperiment"),
	  function(object){
		  assays(object)[[2]]
	  })

## Remove the methods below once the transition from SummarizedExperiment
## to RangedSummarizedExperiment is complete. [H. Pages - May 11, 2015]
setMethod("chromosome", signature(object="SummarizedExperiment"),
	  function(object,...) as.character(seqnames(object)))
setMethod("isSnp", signature(object="SummarizedExperiment"),
	  function(object,...) values(rowRanges(object))$isSnp)
setMethod("lrr", signature(object="SummarizedExperiment"),
	  function(object){
		  assays(object)[[1]]
	  })
setMethod("baf", signature(object="SummarizedExperiment"),
	  function(object){
		  assays(object)[[2]]
	  })

