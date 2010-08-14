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
##SNP-level classes
###########################################################################  
setClass("oligoSnpSet", contains="SnpSet") ## total copy number and genotypes
setClass("CopyNumberSet", contains="eSet") ## total copy number (no genotypes available)

###########################################################################
##Summary-level classes - CNP
###########################################################################
setOldClass("ffdf")
setOldClass("ff_matrix")
setClassUnion("list_or_ffdf", c("list", "ffdf"))

## AssayData elements in AlleleSet are platform dependent.
##
## It is nontrivial to define an initialization method for AlleleSet that can then be extended by
## classes that inherit methods from it.
##
## Easier just to extend SnpSet directly and define accessors for CNSet
##setIs("LinearModelParameter", "AssayData")
##setClass("LinearModelParameter", contains="AssayData")
setClassUnion("LinearModelParameter", c("AssayData", "environment", "list"))
setClass("CNSet", contains = "SnpSuperSet",
	 prototype = prototype(new("VersionedBiobase", versions=c(classVersion("eSet"), CNSet="1.0.0"))))
setClass("CNSetLM", contains="CNSet", representation(lM="list_or_ffdf"))
setMethod("initialize", "CNSetLM", function(.Object, lM=new("list"), ...){
	.Defunct(msg="The CNSetLM class is defunct")	
})
setClass("CNSet", representation(batch="factor", lM="LinearModelParameter"),
	 contains="SnpSet",
	 prototype = prototype(
	 new("VersionedBiobase",
	     versions=c(classVersion("eSet"), CNSet="1.0.1"))))




