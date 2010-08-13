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
setClass("SnpSuperSet", contains=c("AlleleSet", "SnpSet"))

###########################################################################
##Summary-level classes - CNP
###########################################################################
setClass("CNSet", contains="SnpSuperSet")##, representation(batch="factor"))


###########################################################################
##SNP-level classes
###########################################################################  
setClass("oligoSnpSet", contains="SnpSet")
setClass("CopyNumberSet", contains="eSet")

setClassUnion("integerOrMissing", c("integer", "missing", "numeric"))
