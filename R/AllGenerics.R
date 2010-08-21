setGeneric("manufacturer",function(object) standardGeneric("manufacturer"))
setGeneric("manufacturer<-", function(object, value) standardGeneric("manufacturer<-"))
setGeneric("bothStrands", function(object) standardGeneric("bothStrands"))
setGeneric("allele", function(object, allele, strand) standardGeneric("allele"))
setGeneric("annotate", function(object) standardGeneric("annotate"))
setGeneric("getM", function(object) standardGeneric("getM"))
setGeneric("getA", function(object) standardGeneric("getA"))

setGeneric("A", function(object, ...) standardGeneric("A"))
setGeneric("B", function(object, ...) standardGeneric("B"))
setGeneric("A<-", function(object, value) standardGeneric("A<-"))
setGeneric("B<-", function(object, value) standardGeneric("B<-"))

setGeneric("calls<-", function(object, value) standardGeneric("calls<-"))
setGeneric("calls", function(object) standardGeneric("calls"))
setGeneric("confs", function(object, transform=TRUE) standardGeneric("confs"))
setGeneric("confs<-", function(object, value) standardGeneric("confs<-"))
setGeneric("cnConfidence", function(object) standardGeneric("cnConfidence"))
setGeneric("cnConfidence<-", function(object, value) standardGeneric("cnConfidence<-"))
setGeneric("copyNumber", function(object) standardGeneric("copyNumber"))
setGeneric("copyNumber<-", function(object, value) standardGeneric("copyNumber<-"))
setGeneric("chromosome", function(object) standardGeneric("chromosome"))
setGeneric("chromosome<-", function(object, value) standardGeneric("chromosome<-"))
setGeneric("db", function(object) standardGeneric("db"))
setGeneric("kind", function(object) standardGeneric("kind"))
setGeneric("position", function(object) standardGeneric("position"))
setGeneric("isSnp", function(object) standardGeneric("isSnp"))


setGeneric("genomeBuild", function(object) standardGeneric("genomeBuild"))
setGeneric("geometry", function(object) standardGeneric("geometry"))
setGeneric("lM", function(object) standardGeneric("lM"))
setGeneric("lM<-", function(object, value) standardGeneric("lM<-"))
setGeneric("batch", function(object) standardGeneric("batch"))
setGeneric("batch<-", function(object, value) standardGeneric("batch<-"))
setGeneric("batchNames", function(object) standardGeneric("batchNames"))
setGeneric("batchNames<-", function(object,value) standardGeneric("batchNames<-"))


setGeneric("nu", function(object, allele) standardGeneric("nu"))
setGeneric("phi", function(object, allele) standardGeneric("phi"))
setGeneric("sigma2", function(object, allele) standardGeneric("sigma2"))

setGeneric("flags", function(object) standardGeneric("flags"))


setGeneric("batchStatistics", function(object) standardGeneric("batchStatistics"))
setGeneric("batchStatistics<-", function(object,value) standardGeneric("batchStatistics<-"))




