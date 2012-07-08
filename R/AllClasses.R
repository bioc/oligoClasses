###########################################################################
## General DBPDInfo Classes
###########################################################################
setClass("DBPDInfo",
         representation=representation(
           getdb="function",
           tableInfo="data.frame",
           geometry="integer",
           manufacturer="character",
           genomebuild="character",
           annotation="character"))

setClass("SNPPDInfo", contains="DBPDInfo")
setClass("SNPCNVPDInfo", contains="SNPPDInfo")
setClass("ExpressionPDInfo", contains="DBPDInfo")
setClass("TilingPDInfo", contains="DBPDInfo")
setClass("stArrayDBPDInfo", contains="DBPDInfo")
setClass("ExonPDInfo", contains="stArrayDBPDInfo")
setClass("GenePDInfo", contains="stArrayDBPDInfo")

###########################################################################
## Manufacturer-specific PDInfo Classes
###########################################################################
setClass("AffyTilingPDInfo", contains="TilingPDInfo",
         prototype=list(manufacturer="Affymetrix"))
setClass("AffyExpressionPDInfo", contains="ExpressionPDInfo",
         prototype=list(manufacturer="Affymetrix"))

setClass("AffyGenePDInfo", contains="GenePDInfo",
         prototype=list(manufacturer="Affymetrix"))
setClass("AffyExonPDInfo", contains="ExonPDInfo",
         prototype=list(manufacturer="Affymetrix"))
setClass("AffySTPDInfo", contains="AffyExpressionPDInfo")

setClass("AffySNPPDInfo", contains="SNPPDInfo",
         prototype=list(manufacturer="Affymetrix"))
setClass("AffySNPCNVPDInfo", contains="AffySNPPDInfo")

setClass("NgsExpressionPDInfo", contains="ExpressionPDInfo",
         prototype=list(manufacturer="NimbleGen"))
setClass("NgsTilingPDInfo", contains="TilingPDInfo",
         prototype=list(manufacturer="NimbleGen"))

###########################################################################
##Feature-level classes
###########################################################################
setClass("FeatureSet",
         representation=representation(
           manufacturer="character",
           intensityFile="character",
           "VIRTUAL"),
         contains="NChannelSet",
         prototype=prototype(
           manufacturer=NA_character_,
           intensityFile=NA_character_))

setClass("ExpressionFeatureSet", contains="FeatureSet")
setClass("SnpFeatureSet", contains="FeatureSet")
setClass("SnpCnvFeatureSet", contains="SnpFeatureSet")
setClass("TilingFeatureSet", contains="FeatureSet")
setClass("ExonFeatureSet", contains="FeatureSet")
setClass("GeneFeatureSet", contains="FeatureSet")


setClass("AlleleSet", contains="eSet")

###########################################################################
## Combo classes - SNP Summaries - alleles + calls/conf
###########################################################################
## RS is no longer using this class
setClass("SnpSuperSet", contains=c("AlleleSet", "SnpSet"))

###########################################################################
## GenomeAnnotatedDataFrame
###########################################################################
setClass("GenomeAnnotatedDataFrame", contains="AnnotatedDataFrame")

ADF2GDF <- function(object){
	new("GenomeAnnotatedDataFrame",
	    isSnp=as.logical(object$isSnp),
	    position=as.integer(object$position),
	    chromosome=as.integer(object$chromosome),
	    row.names=featureNames(object))
}

setMethod("updateObject", signature(object="GenomeAnnotatedDataFrame"),
	  function(object, ..., verbose=FALSE){
		  ##as(object, "GenomeAnnotatedDataFrame")
		  ADF2GDF(object)
	 })

setMethod("coerce", signature(from="AnnotatedDataFrame", to="GenomeAnnotatedDataFrame"),
	  function(from, to){
		  new("GenomeAnnotatedDataFrame",
		      isSnp=as.logical(from$isSnp),
		      position=as.integer(from$position),
		      chromosome=as.integer(from$chromosome),
		      row.names=featureNames(from))
	  })

###########################################################################
##SNP-level classes
###########################################################################
setClass("gSet", contains="eSet",
	 representation(featureData="GenomeAnnotatedDataFrame",
			genomeBuild="character",
			"VIRTUAL"))

setClass("oligoSnpSet", contains="SnpSet",
	 representation(featureData="GenomeAnnotatedDataFrame")) ## total copy number and genotypes

setClass("CopyNumberSet", contains="gSet") ## total copy number (no genotypes available)
setClass("BeadStudioSet", contains="gSet")


###########################################################################
##Summary-level classes - CNP
###########################################################################
setOldClass("ffdf")
setOldClass("ff_matrix")
setClassUnion("list_or_ffdf", c("list", "ffdf"))
setClassUnion("ff_or_matrix", c("ffdf", "ff_matrix", "matrix"))
## AssayData elements in AlleleSet are platform dependent.
##
## It is nontrivial to define an initialization method for AlleleSet that can then be extended by
## classes that inherit methods from it.
##
## Easier just to extend SnpSet directly and define accessors for CNSet
##setIs("LinearModelParameter", "AssayData")
##setClass("LinearModelParameter", contains="AssayData")
##setClassUnion("LinearModelParameter", c("AssayData", "environment", "list"))
##setClassUnion("NumberGenotype", c("AssayData", "environment", "list"))

setClass("CNSet", contains = "SnpSuperSet",
	 prototype = prototype(new("VersionedBiobase", versions=c(classVersion("eSet"), CNSet="1.0.0"))))

##setClass("GenotypeSummary",
##	 representation(numberGenotype="AssayData",
##			means="AssayData",
##			mads="AssayData"))
####	 prototype=prototype(new("VersionedBiobase", versions=c(GenotypeSummary="1.0.0"))))


setClass("CNSetLM", contains="CNSet", representation(lM="list_or_ffdf"))
setMethod("initialize", "CNSetLM", function(.Object, lM=new("list"), ...){
	.Defunct(msg="The CNSetLM class is defunct")
})

setClass("CNSet", representation(batch="factor",
				 lM="AssayData"),
	 contains="SnpSet",
	 prototype = prototype(
	                       new("VersionedBiobase",
				   versions=c(classVersion("SnpSet"), CNSet="1.0.1"))))
setClass("CNSet", representation(batch="factor",
				 lM="AssayData",
				 numberGenotype="AssayData"),
	 contains="SnpSet",
	 prototype = prototype(
	                       new("VersionedBiobase",
				   versions=c(classVersion("SnpSet"), CNSet="1.0.2"))))
setClass("CNSet", representation(batch="factor",
				 batchStatistics="AssayData"),
	 contains="SnpSet",
	 prototype = prototype(
	                       new("VersionedBiobase",
				   versions=c(classVersion("SnpSet"), CNSet="1.0.3"))))

setClass("CNSet", representation(batch="character",
				 batchStatistics="AssayData"),
	 contains="SnpSet",
	 prototype = prototype(
	                       new("VersionedBiobase",
				   versions=c(classVersion("SnpSet"), CNSet="1.0.4"))))

setClass("CNSet", representation(batch="character",
				 batchStatistics="AssayData",
				 mixtureParams="matrix"),
	 contains="SnpSet",
	 prototype = prototype(
	                       new("VersionedBiobase",
				   versions=c(classVersion("SnpSet"), CNSet="1.0.5"))))

setClass("CNSet", contains="gSet",
	 representation(batch="character",
			batchStatistics="AssayData",
			mixtureParams="ff_or_matrix"),
			##featureData="GenomeAnnotatedDataFrame"),
	 prototype = prototype(
	 new("VersionedBiobase",
	     versions=c(classVersion("SnpSet"), CNSet="1.0.6"))))


setMethod("updateObject", signature(object="CNSet"),
          function(object, ..., verbose=FALSE) {
		  if (verbose) message("updateObject(object = 'CNSet')")
		  obj <- tryCatch(callNextMethod(batch=batch(object)), error=function(e) NULL)
		  if(is.null(obj)){
			  ## must supply batch for batchStatistics to be added
			  if(is(calls(object), "ffdf") | is(calls(object), "ff_matrix"))
				  stopifnot(isPackageLoaded("ff"))
			  if(.hasSlot(object, "mixtureParams")){
				  obj <- new("CNSet",
					     assayData = updateObject(assayData(object),
					     ...., verbose=verbose),
					     phenoData = phenoData(object),
					     experimentData = experimentData(object),
					     annotation = updateObject(annotation(object),
					     ..., verbose=verbose),
					     featureData=updateObject(featureData(object), ..., verbose=verbose),
					     batch=as.character(batch(object)),
					     batchStatistics=batchStatistics(object),
					     mixtureParams=object@mixtureParams)
			  } else {
				  obj <- new("CNSet",
					     assayData = updateObject(assayData(object),
					     ...., verbose=verbose),
					     phenoData = phenoData(object),
					     experimentData = experimentData(object),
					     annotation = updateObject(annotation(object),
					     ..., verbose=verbose),
					     featureData=updateObject(featureData(object), ..., verbose=verbose),
					     batch=as.character(batch(object)),
					     batchStatistics=batchStatistics(object),
					     mixtureParams=matrix(NA, 4, ncol(object)))
			  }
			  if (isCurrent(obj)["CNSet"]) return(obj)
			  return(obj)
		  }
          })


## SetList classes
setClass("eSetList",
	 representation(assayDataList="AssayData",
			phenoData="AnnotatedDataFrame",
			featureDataList="list",
			chromosome="vector",
			annotation="character",
			genomeBuild="character", "VIRTUAL"))
setClass("BeadStudioSetList", contains="eSetList")
setClass("oligoSetList", contains="eSetList")


##---------------------------------------------------------------------------
## classes for ranges
setClass("RangedDataCopyNumber", contains="RangedData",
	 representation("VIRTUAL"))
setClass("RangedDataCNV", contains="RangedDataCopyNumber")
setClass("RangedDataCBS", contains="RangedDataCNV")
setClass("RangedDataHMM", contains="RangedDataCNV")

