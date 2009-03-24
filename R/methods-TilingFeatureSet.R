## setMethod("getM", "TilingQSet", function(object) assayDataElement(object, "M"))

setMethod("getM", "TilingFeatureSet2",
          function(object){
            lc1 <- log2(assayDataElement(object, "channel1"))
            lc2 <- log2(assayDataElement(object, "channel2"))

##             new("TilingFeatureSet2",
##                 logRatio=lc1-lc2,
##                 phenoData=phenoData(object),
##                 experimentData=experimentData(object),
##                 annotation=annotation(object),
##                 featureData=featureData(object))

            lc1-lc2
          })

setMethod("getA", "TilingFeatureSet2",
          function(object){
            lc1 <- log2(assayDataElement(object, "channel1"))
            lc2 <- log2(assayDataElement(object, "channel2"))

##             new("TilingFeatureSet2",
##                 logAverage=(lc1+lc2)/2,
##                 phenoData=phenoData(object),
##                 experimentData=experimentData(object),
##                 annotation=annotation(object),
##                 featureData=featureData(object))

            (lc1+lc2)/2
          })


setMethod("pmPosition", "TilingFeatureSet",
          function(object){
            conn <- db(object)
            tmp <- dbGetQuery(conn, "SELECT fid, position FROM pmfeature")
            tmp <- tmp[order(tmp[["fid"]]),]
            tmp[["position"]]
          })

## setMethod("getX", "TilingFeatureSet",
##           function(object, type){
##             getX(getPD(object), type)
##           })
## 
## setMethod("getY", "TilingFeatureSet",
##           function(object, type){
##             getY(getPD(object), type)
##           })
## 
## setMethod("bgindex", "TilingFeatureSet",
##           function(object){
##             bgindex(getPD(object))
##           })
## 
## setMethod("bg", "TilingFeatureSet",
##           function(object){
##             bgi <- bgindex(object)
##             exprs(object[bgi,])
##           })
