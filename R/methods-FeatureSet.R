## This will contain initialize and accessors to *slots*
## For the moment, every type of computation that depends
##   on annotation, must be in oligo

## Rule of thumb: we should be able to execute the code without further dependencies (only Biobase)

setMethod("initialize", signature(.Object="FeatureSet"), 
          function(.Object,
                   assayData=assayDataNew(exprs=exprs, ...),
                   manufacturer=as.character(NA),
                   platform=as.character(NA),
                   exprs=matrix(numeric(0), nrow=nrow, ncol=ncol),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData=new("MIAME"),
                   annotation=new("character"), ...){
            .Object <- callNextMethod(.Object,
                                      assayData=assayData,
                                      phenoData=phenoData,
                                      experimentData=experimentData,
                                      annotation=annotation,
                                      featureData=featureData)
            .Object@manufacturer <- manufacturer
            .Object@platform <- platform
            .Object
          })


setMethod("platform", signature(object="FeatureSet"),function(object) object@platform)
setReplaceMethod("platform",signature(object="FeatureSet"),
		  function(object,value) {
			  object@platform <- value
			  object
		  })
  
setMethod("manufacturer", signature(object="FeatureSet"), function(object) object@manufacturer)
setReplaceMethod("manufacturer", signature(object="FeatureSet"), 
		function(object, value){
  			object@manufacturer <- value
  			object
})

setMethod("length", signature(x="FeatureSet"),
          function(x) nrow(pData(x)))

setMethod("exprs",
          signature(object="FeatureSet"),
          function(object) assayDataElement(object,"exprs"))
setReplaceMethod("exprs",
		signature(object="FeatureSet"),
		function(object, value) assayDataElementReplace(object, "exprs", value))

setMethod("se.exprs",
          signature(object="FeatureSet"),
          function(object) {
            obj <- assayDataElement(object,"se.exprs")
            if (is.null(obj))
              new("matrix")
            else
              obj
          })
setReplaceMethod("se.exprs",
		signature(object="FeatureSet"),
		function(object, value) assayDataElementReplace(object, "se.exprs", value))
