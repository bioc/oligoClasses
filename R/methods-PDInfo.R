## All methods for DBPDInfo
## Try to use mostly these
##  so extensions will 'just work'
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

setMethod("pmindex", "DBPDInfo",
          function(object, subset=NULL) { 
            if (!is.null(subset))
              warning("Subset not implemented (yet). Returning everything.")
            dbGetQuery(db(object),
                       "select fid from pmfeature")[[1]]
          })

setMethod("mmindex", "DBPDInfo",
          function(object, subset=NULL){
            if (!is.null(subset))
              warning("Subset not implemented (yet). Returning everything.")
            dbGetQuery(db(object),
                       "select fid from mmfeature")[[1]]
          })

setMethod("probeNames", "DBPDInfo",
          function(object, subset=NULL) {
            if (!is.null(subset))
              warning("Subset not implemented (yet). Returning everything.")
            sql <- "select man_fsetid, fid from featureSet, pmfeature where pmfeature.fsetid=featureSet.fsetid"
            tmp <- dbGetQuery(db(object), sql)
            tmp[order(tmp$fid, tmp$man_fsetid), "man_fsetid"]
          })

setMethod("geneNames", "DBPDInfo",
          function(object) {
              unique(probeNames(object))
          })

setMethod("pmSequence", "DBPDInfo",
          function(object){
            sql <- "select seq from sequence, pmfeature where pmfeature.fid=sequence.fid order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })

######## END methods for DBPDInfo

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

setMethod("kind", "ExpressionPDInfo",
          function(object) {
            "expression"
          })

setMethod("kind", "TilingPDInfo",
          function(object) {
            "tiling"
          })


## FIXME: this method should be renamed!
## FIXME: this should query the featureSet table (much faster i think)
setMethod("pmOffset", "AffySNPPDInfo",
          function(object){
            sql <- "select offset from sequence, pmfeature where pmfeature.fid=sequence.fid order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("pmFragmentLength", "AffySNPPDInfo",
          function(object){
            sql <- "select fid, fragment_length from featureSet, pmfeature where pmfeature.fsetid=featureSet.fsetid"
            tmp <- dbGetQuery(db(object), sql)
            idx <- order(tmp[["fid"]])
            tmp[idx, "fragment_length"]
          })

setMethod("pmAllele", "AffySNPPDInfo",
          function(object){
            sql <- "select allele from pmfeature order by fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("pmStrand", "AffySNPPDInfo",
          function(object){
            sql <- "select strand from pmfeature order by fid"
            dbGetQuery(db(object), sql)[[1]]
          })

### For Expression


## FIXME: this method should be renamed!
## FIXME: this should query the featureSet table (much faster i think)
setMethod("pmPosition", "ExpressionPDInfo",
          function(object){
            sql <- "select position from pmfeature order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })


### For Tiling

## FIXME: this method should be renamed!
## FIXME: this should query the featureSet table (much faster i think)
setMethod("pmPosition", "TilingPDInfo",
          function(object){
            sql <- "select position from pmfeature order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })


### AffySNPCNV

## setMethod("pmSequence", "AffySNPCNVPDInfo",
##           function(object, probes.type="snp"){
##             if (probes.type == "snp"){
##               sql <- "select seq from sequence, pmfeature where pmfeature.fid=sequence.fid order by pmfeature.fid"
##             }else if (probes.type == "cnv"){
##               sql <- "SELECT seq FROM sequenceCNV, pmfeatureCNV WHERE pmfeatureCNV.fid=sequenceCNV.fid ORDER BY pmfeatureCNV.fid"
##             }
##             dbGetQuery(db(object), sql)[[1]]
##           })
