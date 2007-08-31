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
                 function(object, value) assayDataElementReplace(object, "cnConfidence", value))


##SnpCallSet
##setMethod("calls", "SnpCallSet", function(object) assayDataElement(object, "calls"))
##setReplaceMethod("calls", signature(object="SnpCallSet", value="matrix"),
##                 function(object, value) assayDataElementReplace(object, "calls", value))
##
##setMethod("callsConfidence", "SnpCallSet", function(object) assayDataElement(object, "callsConfidence"))
##setReplaceMethod("callsConfidence", signature(object="SnpCallSet", value="matrix"),
##                 function(object, value) assayDataElementReplace(object, "callsConfidence", value))

setMethod("db", "SnpCallSet",
          function(object) db(get(annotation(object))))


##SnpCopyNumberSet
##setMethod("copyNumber", "SnpCopyNumberSet", function(object) assayDataElement(object, "copyNumber"))
##setReplaceMethod("copyNumber", signature(object="SnpCopyNumberSet", value="matrix"),
##                 function(object, value) assayDataElementReplace(object, "copyNumber", value))
##
##setMethod("cnConfidence", "SnpCopyNumberSet", function(object) assayDataElement(object, "cnConfidence"))
##setReplaceMethod("cnConfidence", signature(object="SnpCopyNumberSet", value="matrix"),
##                 function(object, value) assayDataElementReplace(object, "cnConfidence", value))


##oligoSnpSet
##setMethod("calls", "oligoSnpSet", function(object) assayDataElement(object, "calls"))
##setReplaceMethod("calls", signature(object="oligoSnpSet", value="matrix"),
##                 function(object, value) assayDataElementReplace(object, "calls", value))
##
##setMethod("callsConfidence", "oligoSnpSet", function(object) assayDataElement(object, "callsConfidence"))
##setReplaceMethod("callsConfidence", signature(object="oligoSnpSet", value="matrix"),
##                 function(object, value) assayDataElementReplace(object, "callsConfidence", value))
##
##
##setMethod("copyNumber", "oligoSnpSet", function(object) assayDataElement(object, "copyNumber"))
##setReplaceMethod("copyNumber", signature(object="oligoSnpSet", value="matrix"),
##                 function(object, value) assayDataElementReplace(object, "copyNumber", value))
##
##setMethod("cnConfidence", "oligoSnpSet", function(object) assayDataElement(object, "cnConfidence"))
##setReplaceMethod("cnConfidence", signature(object="oligoSnpSet", value="matrix"),
##                 function(object, value) assayDataElementReplace(object, "cnConfidence", value))
