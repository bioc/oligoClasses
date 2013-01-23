setGeneric("manufacturer",function(object) standardGeneric("manufacturer"))
setGeneric("manufacturer<-", function(object, value) standardGeneric("manufacturer<-"))
setGeneric("bothStrands", function(object) standardGeneric("bothStrands"))
setGeneric("allele", function(object, allele, strand) standardGeneric("allele"))
setGeneric("annotate", function(object) standardGeneric("annotate"))

## assayData accessors
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
setGeneric("copyNumber", function(object, ...) standardGeneric("copyNumber"))
setGeneric("copyNumber<-", function(object, value) standardGeneric("copyNumber<-"))
setGeneric("baf", function(object) standardGeneric("baf"))
setGeneric("lrr", function(object) standardGeneric("lrr"))
setGeneric("lrr<-", function(object,value) standardGeneric("lrr<-"))
setGeneric("baf<-", function(object,value) standardGeneric("baf<-"))


##GenomeAnnotatedDataFrame accessors
setGeneric("chromosome", function(object, ...) standardGeneric("chromosome"))
setGeneric("chromosome<-", function(object, value) standardGeneric("chromosome<-"))
setGeneric("db", function(object) standardGeneric("db"))
setGeneric("kind", function(object) standardGeneric("kind"))
setGeneric("position", function(object, ...) standardGeneric("position"))
setGeneric("position<-", function(object, value) standardGeneric("position<-"))
setGeneric("isSnp", function(object, ...) standardGeneric("isSnp"))
setGeneric("isSnp<-", function(object, value) standardGeneric("isSnp<-"))
##setGeneric("snpNames", function(object) standardGeneric("snpNames"))
setGeneric("getArm", function(object, ...) standardGeneric("getArm"))

setGeneric("genomeBuild", function(object) standardGeneric("genomeBuild"))
setGeneric("genomeBuild<-", function(object,value) standardGeneric("genomeBuild<-"))
setGeneric("geometry", function(object) standardGeneric("geometry"))


## batchStatistics/copynumber
setGeneric("batch", function(object) standardGeneric("batch"))
setGeneric("batch<-", function(object, value) standardGeneric("batch<-"))
setGeneric("batchNames", function(object) standardGeneric("batchNames"))
setGeneric("batchNames<-", function(object,value) standardGeneric("batchNames<-"))
setGeneric("nu", function(object, allele) standardGeneric("nu"))
setGeneric("phi", function(object, allele) standardGeneric("phi"))
setGeneric("sigma2", function(object, allele) standardGeneric("sigma2"))
setGeneric("batchStatistics", function(object) standardGeneric("batchStatistics"))
setGeneric("batchStatistics<-", function(object,value) standardGeneric("batchStatistics<-"))

setGeneric("flags", function(object) standardGeneric("flags"))
setGeneric("coverage2", function(object) standardGeneric("coverage2"))
setGeneric("state", function(object) standardGeneric("state"))

setGeneric("featuresInRange", function(object, range, FRAME=0, FRAME.LEFT, FRAME.RIGHT, ...)
	   standardGeneric("featuresInRange"))
setGeneric("openff", function(object) standardGeneric("openff"))
setGeneric("closeff", function(object) standardGeneric("closeff"))
setGeneric("checkOrder", function(object, verbose=FALSE) standardGeneric("checkOrder"))
##setGeneric("order2", function(object)  standardGeneric("order2"))
setGeneric("GenomeAnnotatedDataFrameFrom", function(object, annotationPkg, genome="hg19", ...) standardGeneric("GenomeAnnotatedDataFrameFrom"))
setGeneric("numberProbes", function(object) standardGeneric("numberProbes"))
setGeneric("makeFeatureGRanges", function(object, ...) standardGeneric("makeFeatureGRanges"))

## List classes
setGeneric("assayDataList", function(object) standardGeneric("assayDataList"))
setGeneric("featureDataList", function(object) standardGeneric("featureDataList"))

setGeneric("clone", function(object, id, prefix="",...) standardGeneric("clone"))

