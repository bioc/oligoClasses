###########################################################################
##AssayData accessors
###########################################################################
setMethod("calls", "SnpLevelSet", function(object) assayDataElement(object, "calls"))
setReplaceMethod("calls", signature(object="SnpLevelSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "calls", value)) 

setMethod("callsConfidence", "SnpLevelSet", function(object) assayDataElement(object, "callsConfidence"))
setReplaceMethod("callsConfidence", signature(object="SnpLevelSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "callsConfidence", value))

setMethod("copyNumber", "SnpLevelSet", function(object) assayDataElement(object, "copyNumber"))
setReplaceMethod("copyNumber", signature(object="SnpLevelSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "copyNumber", value))

setMethod("cnConfidence", "SnpLevelSet", function(object) assayDataElement(object, "cnConfidence"))
setReplaceMethod("cnConfidence", signature(object="SnpLevelSet", value="matrix"),
                 function(object, value){
                   assayDataElementReplace(object, "cnConfidence", value)
                 })

setMethod("db", "SnpCallSet",
          function(object) db(get(annotation(object))))



