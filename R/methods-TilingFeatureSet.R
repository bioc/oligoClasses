setMethod("getM", "TilingQSet", function(object) assayDataElement(object, "M"))

setMethod("pmPosition", "TilingFeatureSet",
          function(object){
            conn <- db(object)
            tmp <- dbGetQuery(conn, "SELECT fid, position FROM pmfeature")
            tmp <- tmp[order(tmp[["fid"]]),]
            tmp[["position"]]
          })

setMethod("getX", "TilingFeatureSet",
          function(object, type){
            getX(getPD(object), type)
          })

setMethod("getY", "TilingFeatureSet",
          function(object, type){
            getY(getPD(object), type)
          })

setMethod("bgindex", "TilingFeatureSet",
          function(object){
            bgindex(getPD(object))
          })

setMethod("bg", "TilingFeatureSet",
          function(object){
            bgi <- bgindex(object)
            exprs(object[bgi,])
          })
