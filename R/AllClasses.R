###########################################################################
## General DBPDInfo Classes
###########################################################################  
setClass("DBPDInfo",
         representation=representation(
           getdb="function",
           tableInfo="data.frame",
           geometry="integer",
           manufacturer="character",
           genomebuild="character"))

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

setClass("AffyGenePDInfo", contains="GenePDInfo", prototype=list(manufacturer="Affymetrix"))
setClass("AffyExonPDInfo", contains="ExonPDInfo", prototype=list(manufacturer="Affymetrix"))
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
         representation = representation(
           manufacturer="character",
           platform="character",
           "VIRTUAL"),
         contains = "eSet",
         prototype = prototype(
           platform=as.character(NA),
           manufacturer=as.character(NA)))

setClass("ExpressionFeatureSet", contains="FeatureSet")
setClass("SnpFeatureSet", contains="FeatureSet")
setClass("SnpCnvFeatureSet", contains="SnpFeatureSet")
setClass("TilingFeatureSet", contains="FeatureSet")
setClass("ExonFeatureSet", contains="FeatureSet")
setClass("GeneFeatureSet", contains="FeatureSet")


###########################################################################
##Summary-level classes
###########################################################################  
setClass("QuantificationSet", representation("VIRTUAL"), contains="eSet")
setClass("SnpQSet", contains="QuantificationSet")
setClass("SnpCnvQSet", contains="QuantificationSet")
setClass("TilingQSet", contains="QuantificationSet")

###########################################################################
##SNP-level classes
###########################################################################  
setClass("SnpLevelSet", representation("VIRTUAL"), contains="eSet")
setClass("SnpCopyNumberSet", contains = "SnpLevelSet")
setClass("SnpCallSet", contains = "SnpLevelSet")
setClass("oligoSnpSet", contains="SnpLevelSet")

setClass("SnpCallSetPlus", contains=c("SnpQSet",  "SnpCallSet"))
setClass("SnpCnvCallSetPlus", contains=c("SnpCnvQSet", "SnpCallSet"))
