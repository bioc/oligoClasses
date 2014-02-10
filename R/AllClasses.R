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
setClass("HTAPDInfo", contains="stArrayDBPDInfo")

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
setClass("AffyHTAPDInfo", contains="HTAPDInfo")

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
setClass("HTAFeatureSet", contains="FeatureSet")


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

###########################################################################
##SNP-level classes
###########################################################################
setClass("gSet", contains="eSet",
	 representation(##featureData="GenomeAnnotatedDataFrame",
			genome="character",
			"VIRTUAL"))
setClass("SnpSet2", contains="gSet")
##setClass("SnpSet2", contains="SnpSet")
setClass("oligoSnpSet", contains="SnpSet2") ##representation(featureData="GenomeAnnotatedDataFrame"))
setClass("CopyNumberSet", contains="gSet") ## total copy number (no genotypes available)
setClass("BeadStudioSet", contains="gSet")
setClass("BafLrrSet", contains="BeadStudioSet")

#setClass("SomeClass", contains="SnpSet2") ## why will this not work??

###########################################################################
##Summary-level classes - CNP
###########################################################################
setOldClass("ffdf")
setOldClass("ff_matrix")
setClassUnion("list_or_ffdf", c("list", "ffdf"))
setClassUnion("ff_or_matrix", c("ffdf", "ff_matrix", "matrix"))
setClass("CNSet", contains="gSet",
	 representation(batch="character",
			batchStatistics="AssayData",
			mixtureParams="ff_or_matrix",
	                datadir="list"))##,
##	 prototype = prototype(
##	 new("VersionedBiobase",
##	     versions=c(classVersion("SnpSet"), CNSet="1.0.6"))))

setClass("CNSetLM")
setMethod("initialize", "CNSetLM", function(.Object, ...){
	.Defunct(msg="The CNSetLM class is defunct")
})

## SetList classes
setClass("gSetList",
	 representation(assayDataList="AssayData",
			phenoData="AnnotatedDataFrame",
			protocolData="AnnotatedDataFrame",
			experimentData="MIAME",
			featureDataList="list",  ## could be GRangesList...
			chromosome="vector",
			annotation="character",
			genome="character", "VIRTUAL"))
setClass("BeadStudioSetList", contains="gSetList")
setClass("BafLrrSetList", contains="BeadStudioSetList")

setClass("oligoSetList", contains="gSetList")

##---------------------------------------------------------------------------
## classes for ranges
setClass("RangedDataCopyNumber", contains="RangedData",
	 representation("VIRTUAL"))
setClass("RangedDataCNV", contains="RangedDataCopyNumber")
setClass("RangedDataCBS", contains="RangedDataCNV")
setClass("RangedDataHMM", contains="RangedDataCNV")

##setClass("GRangesHMM", contains="GRanges")
##setClass("GRangesHMMList", contains="GRangesList")


##setClass("GRangesList",
##    contains=c("CompressedList", "GenomicRangesList"),
##    representation(
##        unlistData="GRanges",
##        elementMetadata="DataFrame"
##    ),
##    prototype(
##        elementType="GRanges"
##    )
##)

