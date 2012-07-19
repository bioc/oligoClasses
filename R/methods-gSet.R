setMethod("initialize", signature(.Object="gSet"),
          function(.Object,
                   assayData = assayDataNew(...),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   protocolData = phenoData[,integer(0)],
		   genome=c("hg19", "hg18"),
		   featureData,
                   ...) {
		  genome <- match.arg(genome)
		  if(missing(featureData))
			  featureData <- GenomeAnnotatedDataFrameFrom(assayData, annotation, genome=genome)
		  .Object@genome <- genome
		  .Object <- callNextMethod(.Object,
					    assayData = assayData,
					    phenoData = phenoData,
					    featureData = featureData,
					    experimentData = experimentData,
					    annotation = annotation,
					    protocolData = protocolData)
		  return(.Object)
          })

##setValidity("gSet", function(object){
##	if(nrow(object) == 0) return(TRUE)
##	sl <- getSequenceLengths(genomeBuild(object))
##	charChrom <- unique(integer2chromosome(chromosome(object)))
##	fr <- makeFeatureGRanges(object, seqlengths=sl[charChrom])
##	fr2 <-
##})


setMethod("isSnp", signature(object="gSet"),
	  function(object, ...) {
		  isSnp(featureData(object))
	  })

setMethod("isSnp", signature(object="character"),
	  function(object, pkgname, ...){
		  path <- system.file("extdata", package=pkgname)
		  load(file.path(path, "snpProbes.rda"))
		  snpProbes <- get("snpProbes")
		  object %in% snpProbes
	  })

setMethod("db", "gSet",
          function(object) {
		  requireAnnotation(annotation(object)) || stop(paste(annotation(object), "package not available"))
		  get(annotation(object))@getdb()
	  })

setMethod("chromosome", "gSet",
	  function(object, na.rm=FALSE, ...){
		  chromosome(featureData(object), na.rm)
	  })

setReplaceMethod("chromosome", signature(object="gSet", value="integer"),
		 function(object, value){
			 fData(object)$chromosome <-  value
			 object
		 })


setMethod("position", "gSet",
          function(object, na.rm=FALSE, ...){
		  position(featureData(object), na.rm)
          })

setMethod("checkOrder", signature(object="gSet"),
	  function(object, verbose=FALSE){
		  .checkOrder(object, verbose)
	  })


.checkOrder <- function(object, verbose=FALSE){
	d <- diff(order(chromosome(object), position(object)))
	if(any(d < 0)){
		if(verbose)
			warning("Object should be ordered by chromosome and physical position.\n",
				"Try \n",
				"> object <- order(object) \n")
		return(FALSE)
	}
	TRUE
}

chromosomePositionOrder <- function(object, ...){
	is.ordered <- checkOrder(object)
	if(!is.ordered){
		##if(verbose) message("Ordering ", class(object), " object by chromosome and physical position")
		index <- order(chromosome(object), position(object), ...)
		object <- object[index, ]
	}
	return(object)
}

setMethod("genomeBuild", signature(object="gSet"), function(object) object@genome)
setReplaceMethod("genomeBuild", signature(object="gSet", value="character"),
		 function(object, value){
			 object@genome <- value
			 return(object)
		 })

setMethod("show", signature(object="gSet"),
	  function(object){
		  callNextMethod(object)
		  cat("genome: ", genomeBuild(object), "\n")
	  })


setMethod("makeFeatureGRanges", signature(object="gSet"),
	  function(object, ...){
		  makeFeatureGRanges(featureData(object), genomeBuild(object))
##		  sl <- getSequenceLengths(genomeBuild(object))
##		  gr <- GRanges(paste("chr", chromosome(object), sep=""),
##				IRanges(position(object), width=1))
##		  seqlengths(gr) <- sl[match(unique(seqnames(gr)), names(sl))]
##		  return(gr)
	  })

setMethod("getArm", signature(object="gSet"), function(object, ...){
	.getArm(chromosome(object), position(object), genomeBuild(object))
})
