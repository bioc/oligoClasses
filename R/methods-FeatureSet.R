## This will contain initialize and accessors to *slots*
## For the moment, every type of computation that depends
##   on annotation, must be in oligo

## Rule of thumb: we should be able to execute the code without further dependencies (only Biobase)

setMethod("initialize", signature(.Object="FeatureSet"), 
          function(.Object,
                   assayData=assayDataNew(...),
                   manufacturer=as.character(NA),
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
            .Object
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

setMethod("kind",
          signature(object="FeatureSet"),
          function(object){
            kind(getPlatformDesign(object))
          })

setMethod("getPlatformDesign",
          signature(object= "FeatureSet"),
          function(object){
            pdn <- annotation(object)
            library(pdn,character.only=TRUE)
            return(get(pdn,pos=paste("package:",pdn,sep="")))
          })
getPD <- getPlatformDesign

setMethod("bgSequence",
          signature(object="FeatureSet"),
          function(object){
            bgSequence(getPlatformDesign(object))
          })

setMethod("pmSequence",
          signature(object="FeatureSet"),
          function(object){
            pmSequence(getPlatformDesign(object))
          })

setMethod("mmSequence",
          signature(object="FeatureSet"),
          function(object){
            mmSequence(getPlatformDesign(object))
          })

setMethod("genomeBuild",
          signature(object="FeatureSet"),
          function(object){
            genomeBuild(getPlatformDesign(object))
          })

setMethod("pmChr", "FeatureSet",
          function(object){
            conn <- db(object)
            if (is(object, "TilingFeatureSet") & manufacturer(object) == "Affymetrix"){
              sql <- paste("SELECT fid, chrom_id as chrom",
                           "FROM pmfeature",
                           "INNER JOIN chrom_dict",
                           "USING(chrom)")
            }else{
              sql <- paste("SELECT fid, chrom",
                           "FROM pmfeature, featureSet",
                           "WHERE pmfeature.fsetid=featureSet.fsetid")
            }
            tmp <- dbGetQuery(conn, sql)
            tmp <- tmp[order(tmp[["fid"]]),]
            tmp[["chrom"]]
          })

setMethod("db", "FeatureSet",
          function(object){
            db(getPlatformDesign(object))
          })

setMethod("geneNames", "FeatureSet",
          function(object){
            geneNames(getPlatformDesign(object))
          })

setMethod("pmindex", "FeatureSet",
          function(object, subset=NULL){
            pmindex(getPlatformDesign(object), subset=subset)
          })

setMethod("mmindex", "FeatureSet",
          function(object, subset=NULL){
            mmindex(getPD(object), subset=subset)
          })

setMethod("probeNames", "FeatureSet",
          function(object, subset=NULL) {
            if (!is.null(subset))
              warning("ignoring subset arg, feature not implemented")
            probeNames(getPlatformDesign(object))
          })

setMethod("bgindex", "FeatureSet",
          function(object, subset=NULL){
            bgindex(getPD(object), subset=subset)
          })

setMethod("pmPosition", "FeatureSet",
          function(object){
            conn <- db(object)
            sql <- paste("SELECT fid, position",
                         "FROM pmfeature",
                         "INNER JOIN featureSet USING(fsetid)")
            tmp <- dbGetQuery(conn, sql)
            tmp <- tmp[order(tmp[["fid"]]),]
            tmp[["position"]]
          })
