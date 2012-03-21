## Only accessors for FeatureSet objects

## In addition to the slots inherited from eSet, FeatureSet has:
##   - manufacturer
##   - 'exprs' in assayData
##   - 'channel1' and 'channel2' in assayData

setMethod("manufacturer", signature(object="FeatureSet"),
          function(object) object@manufacturer)
setReplaceMethod("manufacturer", signature(object="FeatureSet"),
		function(object, value){
                  object@manufacturer <- value
                  object
})

setMethod("exprs",
          signature(object="FeatureSet"),
          function(object)
          assayDataElement(object,"exprs"))
setReplaceMethod("exprs",
		signature(object="FeatureSet"),
		function(object, value)
                 assayDataElementReplace(object, "exprs", value))

setMethod("db", "FeatureSet",
          function(object){
            db(get(annotation(object)))
          })

setMethod("kind", "FeatureSet",
          function(object){
            kind(get(annotation(object)))
          })

setMethod("geometry", "FeatureSet",
          function(object)
          geometry(getPD(object))
          )
