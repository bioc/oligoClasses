setMethod("initialize", signature(.Object="gSetList"),
	  function(.Object,
		   assayDataList=AssayDataList(...),
		   featureDataList=GenomeAnnotatedDataFrameFromList(assayDataList),
		   chromosome=integer(),
		   phenoData=annotatedDataFrameFrom(assayDataList, byrow=FALSE),
		   protocolData=phenoData[, integer(0)],
		   experimentData=new("MIAME"),
		   annotation=character(),
		   genome=character(),
		   ...){
		  callNextMethod(.Object,
				 assayDataList=assayDataList,
				 featureDataList=featureDataList,
				 phenoData=phenoData,
				 chromosome=chromosome,
				 annotation=annotation,
				 genome=genome,
				 protocolData=protocolData,
				 experimentData=experimentData,
				 ...)
	  })
setMethod("annotation", signature(object="gSetList"), function(object) object@annotation)
setMethod("genomeBuild", signature(object="gSetList"), function(object) object@genome)
setReplaceMethod("genomeBuild", signature(object="gSetList", value="character"), function(object, value){
	object@genome <- value
	object
})
