setMethod("initialize", "CopyNumberSet",
          function(.Object,
                   assayData = assayDataNew(copyNumber = copyNumber,
                                            cnConfidence = cnConfidence, ...),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   protocolData = phenoData[,integer(0)],
                   copyNumber = new("matrix"),
                   cnConfidence = matrix(numeric(),
                                            nrow=nrow(copyNumber), ncol=ncol(copyNumber),
                                            dimnames=dimnames(copyNumber)),
                   ...) {
		  .Object <- callNextMethod(.Object,
					    assayData = assayData,
					    phenoData = phenoData,
					    featureData = featureData,
					    experimentData = experimentData,
					    annotation = annotation,
					    protocolData = protocolData)
		  if(checkAnnotation(annotation))
			  .Object <- annotate(.Object)
		  return(.Object)
          })

setMethod("initialize", "oligoSnpSet",
	  function(.Object,
		   call=new("matrix"),
		   callProbability=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
		   copyNumber=matrix(numeric(), nrow=nrow(call), ncol=ncol(call),  dimnames=dimnames(call)),
		   cnConfidence=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
		   annotation=character(), ...){
		  .Object <- callNextMethod(.Object,
					    call=call,
					    callProbability=callProbability,
					    copyNumber=copyNumber,
					    cnConfidence=cnConfidence,
					    annotation=annotation, ... )
		  if(checkAnnotation(annotation))
			  .Object <- annotate(.Object)  		  
		  return(.Object)
	  })

setValidity("oligoSnpSet", function(object) {
	assayDataValidMembers(assayData(object), c("call", "callProbability", "copyNumber", "cnConfidence"))
})

## RS: ask BC about this... initialization method for CNSet does not work when this is uncommented
##setValidity("AlleleSet",
##            function(object){
##              grp1 <- c("alleleA", "alleleB")
##              grp2 <- c("senseAlleleA", "senseAlleleB",
##                        "antisenseAlleleA", "antisenseAlleleB")
##              elem <- assayDataElementNames(object)
##              ok <- all(grp1 %in% elem) || all(grp2 %in% elem)
##              f <- function(x) paste("'", x, "'", collapse=" + ", sep="")
##              if (!ok){
##                paste("Elements of 'AlleleSummarySet' must be:",
##                      f(grp1), "OR", f(grp2))
##              }else{
##                TRUE
##              }
##            })

setMethod("initialize", "CNSet",
          function(.Object,
		   call=new("matrix"),		   		   
		   CA=matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
		   CB=matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
                   callProbability=matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
                   alleleA = matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
                   alleleB = matrix(NA, nrow(call), ncol(call), dimnames=dimnames(call)),
                   phenoData = annotatedDataFrameFrom(call, byrow=FALSE),		   
		   featureData=annotatedDataFrameFrom(call, byrow=TRUE),
		   experimentData=new("MIAME"),
		   protocolData=phenoData[, integer(0)],
		   annotation=character(), ... ){
		  ##The ... should be additional assayDataElements (e.g., for a class that extends CNSet)
		  .Object <- callNextMethod(.Object,
					    call=call,
					    callProbability=callProbability,
					    alleleA=alleleA,
					    alleleB=alleleB,
					    CA=CA,
					    CB=CB,
					    phenoData=phenoData,
					    featureData=featureData,
					    experimentData=experimentData,
					    protocolData=protocolData,
					    annotation=annotation,
					    ...)
		  if(checkAnnotation(annotation))
			  .Object <- annotate(.Object)
		  return(.Object)
          })


setValidity("CNSet", function(object) {
	assayDataValidMembers(assayData(object), c("CA", "CB", "call", "callProbability", "alleleA", "alleleB"))
})
