setMethod("annotation", "DBPDInfo",
          function(object) object@annotation)

setMethod("initialize", "DBPDInfo",
          function(.Object, ...) {
            .Object <- callNextMethod()
            tInfo <- dbGetQuery(db(.Object), "select * from table_info")
            .Object@tableInfo <- tInfo
            .Object
          })

setMethod("manufacturer", "DBPDInfo",
          function(object) object@manufacturer)

setMethod("genomeBuild", "DBPDInfo",
          function(object) object@genomebuild)

setMethod("geometry", "DBPDInfo",
          function(object) object@geometry)

setMethod("db", signature(object="DBPDInfo"),
          function(object) object@getdb())

setMethod("kind", "AffySNPPDInfo",
          function(object) {
            "SNP"
          })
  
setMethod("kind", "AffyExpressionPDInfo",
          function(object) {
              "expression"
          })

setMethod("kind", "AffySNPCNVPDInfo",
          function(object) {
              "SNPCNV"
          })
 
setMethod("kind", "AffyGenePDInfo",
          function(object) {
              "gene"
          })

setMethod("kind", "AffyExonPDInfo",
          function(object) {
            "exon"
          })

setMethod("kind", "AffyHTAPDInfo",
          function(object) {
            "hta"
          })

setMethod("kind", "ExpressionPDInfo",
          function(object) {
            "expression"
          })

setMethod("kind", "TilingPDInfo",
          function(object) {
            "tiling"
          })
