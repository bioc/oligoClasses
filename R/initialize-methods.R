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
	lm <- assayDataNew(tau2A=initializeBigMatrix("tau2A", nr, nc),
			   tau2B=initializeBigMatrix("tau2B", nr, nc),
			   sig2A=initializeBigMatrix("sig2A", nr, nc),
			   sig2B=initializeBigMatrix("sig2B", nr, nc),
			   nuA=initializeBigMatrix("nuA", nr, nc),
			   nuB=initializeBigMatrix("nuB", nr, nc),
			   phiA=initializeBigMatrix("phiA", nr, nc),
			   phiB=initializeBigMatrix("phiB", nr, nc),
			   phiPrimeA=initializeBigMatrix("phiPrimeA", nr, nc),
			   phiPrimeB=initializeBigMatrix("phiPrimeB", nr, nc),
			   corrAB=initializeBigMatrix("corrAB", nr, nc),
			   corrBB=initializeBigMatrix("corrBB", nr, nc),
			   corrAA=initializeBigMatrix("corrAA", nr, nc),
			   flags=initializeBigMatrix("flags", nr, nc))
	return(lm)
}

initializeNumberGenotypeFrom <- function(object, batchnames, nr){
	if(missing(nr)) nr <- nrow(object)
	nc <- 3
	nGt <- vector("list", length(batchnames))
	elem.names <- paste("N_", batchnames, sep="")
	for(i in seq_along(batchnames)) nGt[[i]] <- initializeBigMatrix(elem.names[i], nr, nc)
	names(nGt) <- batchnames
	nGt <- lapply(nGt, function(x){ colnames(x) <- c("AA", "AB", "BB"); return(x)})
	aD <- do.call(assayDataNew, nGt)
	return(aD)
}


setMethod("initialize", signature=signature(.Object="GenotypeSummary"),
	  function(.Object, numberGenotypes, means, mads){
		  .Object@numberGenotypes <- numberGenotypes
		  .Object@means <- means
		  .Object@mads <- mads
	  })

setMethod("initialize", "CNSet",
	  function(.Object, lM, batch, numberGenotype, ...){
		  .Object@numberGenotype <- assayDataNew()
		  .Object@lM <- assayDataNew()
		  .Object <- callNextMethod(.Object, ...)
		  if(missing(batch)){
			  stop("Must specify factor 'batch'. See ?CNSet-class for details.")
		  } else .Object@batch <- batch
		  if(missing(lM)){
			  lM(.Object) <- initializeLmFrom(.Object)
		  } else lM(.Object) <- lM
		  batchNames(.Object) <- unique(as.character(batch))
		  if(missing(numberGenotype)){
			  .Object@numberGenotype <- initializeNumberGenotypeFrom(.Object, batchNames(.Object))
		  } else .Object@numberGenotype <- numberGenotype
##		  if(missing(genotypeSummary)){
##			  tmp <- initializeGenotypeSummaryFrom(.Object)
##			  gtSum <- new("GenotypeSummary",
##				       numberGenotypes=tmp[["numberGenotypes"]],
##				       means=tmp[["means"]],
##				       mads=tmp[["mads"]])
##			  .Object@genotypeSummary <- gtSum
##		  } else .Object@genotypeSummary <- genotypeSummary
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





