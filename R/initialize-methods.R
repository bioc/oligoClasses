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


setAs("CNSet", "CopyNumberSet",
      function(from){
	      new("CopyNumberSet",
		  copyNumber=totalCopynumber(from, i=1:nrow(from), j=1:ncol(from)),
		  annotation=annotation(from),
		  featureData=featureData(from),
		  phenoData=phenoData(from),
		  experimentData=experimentData(from),
		  protocolData=protocolData(from))
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
setValidity("AlleleSet",
            function(object){
              grp1 <- c("alleleA", "alleleB")
              grp2 <- c("senseAlleleA", "senseAlleleB",
                        "antisenseAlleleA", "antisenseAlleleB")
              elem <- assayDataElementNames(object)
              ok <- all(grp1 %in% elem) || all(grp2 %in% elem)
              f <- function(x) paste("'", x, "'", collapse=" + ", sep="")
              if (!ok){
                paste("Elements of 'AlleleSet' must be:",
                      f(grp1), "OR", f(grp2))
              }else{
                TRUE
              }
            })
setMethod("initialize", "SnpSuperSet", function(.Object,  ...) callNextMethod(.Object, ...))


initializeLmFrom <- function(object){
	nr <- nrow(object)
	nc <- length(unique(batch(object)))
	lm <- assayDataNew(N.AA=initializeBigMatrix("N.AA", nr, nc),
			   N.AB=initializeBigMatrix("N.AB", nr, nc),
			   N.BB=initializeBigMatrix("N.BB", nr, nc),
			   medianA.AA=initializeBigMatrix("medianA.AA", nr, nc),
			   medianA.AB=initializeBigMatrix("medianA.AB", nr, nc),
			   medianA.BB=initializeBigMatrix("medianA.BB", nr, nc),
			   medianB.AA=initializeBigMatrix("medianB.AA", nr, nc),
			   medianB.AB=initializeBigMatrix("medianB.AB", nr, nc),
			   medianB.BB=initializeBigMatrix("medianB.BB", nr, nc),
			   madA.AA=initializeBigMatrix("madA.AA", nr, nc, vmode="double"),
			   madA.AB=initializeBigMatrix("madA.AB", nr, nc, vmode="double"),
			   madA.BB=initializeBigMatrix("madA.BB", nr, nc, vmode="double"),
			   madB.AA=initializeBigMatrix("madB.AA", nr, nc, vmode="double"),
			   madB.AB=initializeBigMatrix("madB.AB", nr, nc, vmode="double"),
			   madB.BB=initializeBigMatrix("madB.BB", nr, nc, vmode="double"),
			   tau2A.AA=initializeBigMatrix("tau2A.AA", nr, nc, vmode="double"),
			   tau2A.BB=initializeBigMatrix("tau2A.BB", nr, nc, vmode="double"),
			   tau2B.AA=initializeBigMatrix("tau2B.AA", nr, nc, vmode="double"),
			   tau2B.BB=initializeBigMatrix("tau2B.BB", nr, nc, vmode="double"),
			   nuA=initializeBigMatrix("nuA", nr, nc, vmode="double"),
			   nuB=initializeBigMatrix("nuB", nr, nc, vmode="double"),
			   phiA=initializeBigMatrix("phiA", nr, nc, vmode="double"),
			   phiB=initializeBigMatrix("phiB", nr, nc, vmode="double"),
			   phiPrimeA=initializeBigMatrix("phiPrimeA", nr, nc, vmode="double"),
			   phiPrimeB=initializeBigMatrix("phiPrimeB", nr, nc, vmode="double"),
			   corrAB=initializeBigMatrix("corrAB", nr, nc, vmode="double"),
			   corrBB=initializeBigMatrix("corrBB", nr, nc, vmode="double"),
			   corrAA=initializeBigMatrix("corrAA", nr, nc, vmode="double"),
			   flags=initializeBigMatrix("flags", nr, nc))
	return(lm)
}


setMethod("initialize", "CNSet",
	  function(.Object, batchStatistics, batch, ...){
##		  if(missing(featureData) & missing(annotationPackage))
##			  stop("must specify valid annotation package.  See annotationPackages()")
##		  if(missing(featureData)) {
##			  featureData <- featureDataFrom(annotationPackage)
##		  } else stopifnot(all(c("chromosome", "isSnp", "position") %in% varLabels(featureData)))
		  .Object@batchStatistics <- assayDataNew()
		  .Object <- callNextMethod(.Object, ...)
##		  .Object@featureData <- featureData
		  if(missing(batch)){
			  stop("Must specify factor 'batch'. See ?CNSet-class for details.")
		  } else .Object@batch <- batch
		  if(missing(batchStatistics)){
			  batchStatistics(.Object) <- initializeLmFrom(.Object)
		  } else batchStatistics(.Object) <- batchStatistics
		  batchNames(.Object) <- unique(as.character(batch))
		  return(.Object)
})

setValidity("CNSet", function(object){
	if(!assayDataValidMembers(assayData(object), c("alleleA", "alleleB", "call", "callProbability"))){
		message("assay data members must be 'alleleA', 'alleleB', 'call', 'callProbability'")
		return(FALSE)
	}
	if(length(batch) == ncol(object)){
		message("Factor 'batch' must be the same length as the number of samples.  See ?CNSet-class for details")
		return(FALSE)
	}
})

initializeGenotypeSummaryFrom <- function(object){
	nr <- nrow(object)
	nc <- 3
	bns <- batchNames(object)
	elem.names <- paste("N_", bns, sep="")
	nGt <- vector("list", length(bns))
	for(i in seq_along(bns)) nGt[[i]] <- initializeBigMatrix(elem.names[i], nr, nc)
	names(nGt) <- elem.names
	numberGt <- do.call(assayDataNew, nGt)

	elem.names <- paste("mns_", bns, sep="")
	mns <- vector("list", length(bns))
	for(i in seq_along(bns)) mns[[i]] <- initializeBigMatrix(elem.names[i], nr, nc)
	mns <- do.call(assayDataNew, mns)

	elem.names <- paste("mads_", bns, sep="")
	mads <- vector("list", length(bns))
	for(i in seq_along(bns)) mads[[i]] <- initializeBigMatrix(elem.names[i], nr, nc)
	mads <- do.call(assayDataNew, mads)
	return(list(numberGenotypes=numberGt, means=mns, mads=mads))
}





