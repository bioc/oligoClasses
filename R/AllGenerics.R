setGeneric("platform", function(object) standardGeneric("platform"))
setGeneric("platform<-", function(object, value) standardGeneric("platform<-"))
setGeneric("manufacturer",function(object) standardGeneric("manufacturer"))
setGeneric("manufacturer<-",
                      function(object, value) standardGeneric("manufacturer<-"))
setGeneric("senseThetaA", function(object) standardGeneric("senseThetaA"))
setGeneric("senseThetaB", function(object) standardGeneric("senseThetaB"))
setGeneric("antisenseThetaA", function(object) standardGeneric("antisenseThetaA"))
setGeneric("antisenseThetaB", function(object) standardGeneric("antisenseThetaB"))
setGeneric("getM", function(object) standardGeneric("getM"))
setGeneric("getA", function(object) standardGeneric("getA"))
setGeneric("plotM", function(object, i, ...) standardGeneric("plotM"))
setGeneric("thetaA", function(object) standardGeneric("thetaA"))
setGeneric("thetaB", function(object) standardGeneric("thetaB"))


setGeneric("calculateCopyNumber", function(object, ...) standardGeneric("calculateCopyNumber"))
setGeneric("calls<-", function(object, value) standardGeneric("calls<-"))
setGeneric("callsConfidence<-",
           function(object, value) standardGeneric("callsConfidence<-"))
setGeneric("calls", function(object) standardGeneric("calls"))
setGeneric("callsConfidence", function(object) standardGeneric("callsConfidence"))
setGeneric("copyNumber<-",
           function(object, value) standardGeneric("copyNumber<-"))
setGeneric("cnConfidence<-",
           function(object, value) standardGeneric("cnConfidence<-"))
setGeneric("copyNumber", function(object) standardGeneric("copyNumber"))
setGeneric("cnConfidence", function(object) standardGeneric("cnConfidence"))
setGeneric("chromosome", function(object) standardGeneric("chromosome"))
setGeneric("chromosome<-", function(object, value) standardGeneric("chromosome<-"))
setGeneric("db", function(object) standardGeneric("db"))
setGeneric("position", function(object) standardGeneric("position"))
