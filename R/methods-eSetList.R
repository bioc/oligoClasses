setMethod("initialize", signature(.Object="eSetList"),
	  function(.Object,
		   assayDataList=AssayDataList(...),
		   featureDataList=GenomeAnnotatedDataFrameFromList(assayDataList),
		   chromosome=integer(),
		   phenoData=annotatedDataFrameFrom(assayDataList, byrow=FALSE),
		   annotation=character(),
		   genome=character(),
		   ...){
		  callNextMethod(.Object,
				 assayDataList=assayDataList,
				 featureDataList=featureDataList,
				 phenoData=phenoData,
				 chromosome=chromosome,
				 annotation=annotation,
				 genome=genome, ...)
	  })
setMethod("annotation", signature(object="eSetList"), function(object) object@annotation)
setMethod("genomeBuild", signature(object="eSetList"), function(object) object@genome)
