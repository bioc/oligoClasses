setMethod("initialize", signature(.Object="eSetList"),
	  function(.Object,
		   assayDataList=AssayDataList(...),
		   featureDataList=GenomeAnnotatedDataFrameFromList(assayDataList),
		   chromosome=integer(),
		   phenoData=annotatedDataFrameFrom(assayDataList, byrow=FALSE),
		   annotation=character(),
		   genomeBuild=character(),
		   ...){
		  callNextMethod(.Object,
				 assayDataList=assayDataList,
				 featureDataList=featureDataList,
				 phenoData=phenoData,
				 chromosome=chromosome,
				 annotation=annotation,
				 genomeBuild=genomeBuild)
	  })
