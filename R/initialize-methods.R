setMethod("initialize", signature(.Object="CopyNumberSet"),
          function(.Object,
                   assayData = assayDataNew(copyNumber = copyNumber,
                                            cnConfidence = cnConfidence, ...),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   protocolData = phenoData[,integer(0)],
                   copyNumber = new("matrix"),
                   cnConfidence = matrix(numeric(),
                                            nrow=nrow(copyNumber), ncol=ncol(copyNumber),
                                            dimnames=dimnames(copyNumber)),
		   featureData=GenomeAnnotatedDataFrameFrom(assayData, annotation, genome=genome),
		   genome=c("hg19", "hg18"),
                   ...) {
		  if(nrow(copyNumber)>0){
			  if(!is(copyNumber[, 1], "integer")) stop("copyNumber should be supplied as a matrix of integers (original scale * 100). See integerMatrix in the oligoClasses package for the conversion to integer matrices")
		  }
		  .Object <- callNextMethod(.Object,
					    assayData = assayData,
					    phenoData = phenoData,
					    featureData = featureData,
					    experimentData = experimentData,
					    annotation = annotation,
					    protocolData = protocolData)
		  if(nrow(assayData[["copyNumber"]]) == 0){
			  .Object@genome <- ""
		  } else{
			  .Object@genome <- match.arg(genome)
		  }
		  return(.Object)
          })

setValidity("CopyNumberSet", function(object){
	assayDataValidMembers(assayData(object), c("copyNumber", "cnConfidence"))
	msg <- isValidGenomeAnnotatedDataFrame(featureData(object))
	if(is.character(msg)) return(msg) else TRUE
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

##setMethod("initialize", "oligoSnpSet",
##	  function(.Object,
##		   call=new("matrix"),
##		   callProbability=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
##		   copyNumber=matrix(numeric(), nrow=nrow(call), ncol=ncol(call),  dimnames=dimnames(call)),
##		   ##cnConfidence=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
##		   assayData=Biobase::assayDataNew(call=call,
##		   callProbability=callProbability,
##		   copyNumber=copyNumber, ...),
##		   annotation=character(),
##		   phenoData,
##		   featureData, ##=GenomeAnnotatedDataFrameFrom(call, annotation),
##		   experimentData,
##		   protocolData,
##		   genome=c("hg19", "hg18"),
##		   ...){
##		  if(nrow(copyNumber)>0){
##			  if(!is(copyNumber[, 1], "integer")) stop("copyNumber should be supplied as a matrix of integers (original scale * 100). See integerMatrix in the oligoClasses package for the conversion to integer matrices")
##		  }
##		  nms <- names(list(...))
##		  if(length(nms) > 0){
##			  ## check that each element in ... is a matrix of integers
##			  for(i in seq_along(nms)){
##				  elt <- list(...)[[nms[i]]]
##				  if(nrow(elt) > 0)
##					  if(!is(elt[,1], "integer")) stop("all assay data elements must be integers. For copy number, use original scale * 100 and for B allele frequencies use original scale * 1000. See integerMatrix in the oligoClasses package for the conversion to integer matrices")
##			  }
##		  }
##		  if(missing(featureData))
##			  featureData <- GenomeAnnotatedDataFrameFrom(assayData, annotation)
##		  if(missing(phenoData))
##			  phenoData <- Biobase::annotatedDataFrameFrom(call, byrow=FALSE)
##		  if(missing(experimentData))
##			  experimentData <- new("MIAME")
##		  if(missing(protocolData))
##			  protocolData <- phenoData[, integer(0)]
##		  .Object@genome <- match.arg(genome)
##		  .Object <- callNextMethod(.Object,
##					    assayData=assayData,
##					    annotation=annotation,
##					    featureData=featureData,
##					    experimentData=experimentData,
##					    phenoData=phenoData,
##					    protocolData=protocolData,
##					    ...)
##		  return(.Object)
##	  })

##harmomonizeAssayData <- function(assayData){
##	nms <- names(ls(assayData))
##	nr <- rep(NA, length(nms))
##	for(i in seq_along(nms)){
##		elt <- assayData[[nms[i]]]
##		nr[i] <- nrow(elt)
##	}
##}

setMethod("initialize", "oligoSnpSet",
	  function(.Object,
		   call=new("matrix"),
		   callProbability=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
		   copyNumber=matrix(numeric(), nrow=nrow(call), ncol=ncol(call),  dimnames=dimnames(call)),
		   ##cnConfidence=matrix(numeric(), nrow=nrow(call), ncol=ncol(call), dimnames=dimnames(call)),
		   assayData=Biobase::assayDataNew(call=call,
		   callProbability=callProbability,
		   copyNumber=copyNumber, ...),
		   annotation=character(),
		   phenoData,
		   featureData, ##=GenomeAnnotatedDataFrameFrom(call, annotation),
		   experimentData,
		   protocolData,
		   genome=c("hg19", "hg18"),
		   ...){
		  if(nrow(copyNumber)>0){
			  if(!is(copyNumber[, 1], "integer")) stop("copyNumber should be supplied as a matrix of integers (original scale * 100). See integerMatrix in the oligoClasses package for the conversion to integer matrices")
		  }
		  nms <- names(list(...))
		  if(length(nms) > 0){
			  ## check that each element in ... is a matrix of integers
			  for(i in seq_along(nms)){
				  elt <- list(...)[[nms[i]]]
				  if(nrow(elt) > 0)
					  if(!is(elt[,1], "integer")) stop("all assay data elements must be integers. For copy number, use original scale * 100 and for B allele frequencies use original scale * 1000. See integerMatrix in the oligoClasses package for the conversion to integer matrices")
			  }
		  }
		  genome <- match.arg(genome)
		  if(missing(featureData))
			  featureData <- GenomeAnnotatedDataFrameFrom(assayData, annotation, genome=genome)
		  if(missing(phenoData))
			  phenoData <- Biobase::annotatedDataFrameFrom(call, byrow=FALSE)
		  if(missing(experimentData))
			  experimentData <- new("MIAME")
		  if(missing(protocolData))
			  protocolData <- phenoData[, integer(0)]
		  .Object <- callNextMethod(.Object,
					    assayData=assayData,
					    annotation=annotation,
					    featureData=featureData,
					    experimentData=experimentData,
					    phenoData=phenoData,
					    protocolData=protocolData,
					    ...)
		  .Object@genome <- genome
		  return(.Object)
	  })

setValidity("oligoSnpSet", function(object){
	##nms <- ls(assayData(object))
	Biobase::assayDataValidMembers(assayData(object), c("call", "callProbability", "copyNumber"))
	msg <- isValidGenomeAnnotatedDataFrame(featureData(object))
	if(nrow(copyNumber(object)) > 0){
		if(!is.integer(copyNumber(object)[,1])) return("copyNumber should be a matrix of integers (original scale * 100). Use integerMatrix(x, 100) for converting 'x' to a matrix of integers.")
	}
	if("baf" %in% ls(assayData(object))){
		b <- assayData(object)[["baf"]]
		if(nrow(b) > 0){
			if(!is.integer(b[,1])) return("B allele frequencies should be a matrix of integers (original scale * 1000). See integerMatrix(x, 1000) for converting 'x' to a matrix of integers.")
		}
	}
	if(nrow(object) > 0){
		genome <- genomeBuild(object)
		if(!genome %in% c("hg18", "hg19")) return("Supported values for genome are hg18 and hg19")
	}
	if(is.character(msg)) return(msg)
	validObject(phenoData(object))
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
	if(nc > 1) nc <- nc+1 ## add extra column for grand mean
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

initializeLmFrom2 <- function(object, batch){
	nr <- nrow(object)
	nc <- length(unique(batch))
	if(nc > 1) nc <- nc+1 ## add extra column for grand mean
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
	if(nc > 1) {
		sampleNames(lm) <- c(unique(batch), "grandMean")
	} else sampleNames(lm) <- unique(batch)
	return(lm)
}


setMethod("initialize", "CNSet",
	  function(.Object,
		   alleleA=new("matrix"),
		   alleleB=alleleA,
		   call=alleleA,
		   callProbability=alleleA,
		   assayData=assayDataNew(alleleA=alleleA,
		                          alleleB=alleleB,
		                          call=call,
		                          callProbability=callProbability, ...),
		   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
		   protocolData=phenoData[, integer(0)],
		   experimentData=new("MIAME"),
		   annotation=character(),
		   featureData,
		   batch=character(ncol(alleleA)),
		   batchStatistics=initializeLmFrom2(alleleA, batch),
		   genome=c("hg19", "hg18"),
		   mixtureParams=new("matrix"),
		   datadir=list("", integer(), ""), ...){
		  genome <- match.arg(genome)
		  if(missing(featureData))
			  featureData <- GenomeAnnotatedDataFrameFrom(assayData, annotation, genome=genome)
		  .Object@mixtureParams <- mixtureParams
		  .Object@batch <- batch
		  .Object@batchStatistics <- batchStatistics
		  .Object <- callNextMethod(.Object,
					    assayData=assayData,
					    phenoData=phenoData,
					    featureData=featureData,
					    experimentData=experimentData,
					    annotation=annotation,
					    protocolData=protocolData, ...)
		  .Object@datadir <- datadir
		  .Object@genome <- genome
		  if(nrow(.Object)==0) .Object@genome <- character()
		  return(.Object)
})

setValidity("CNSet", function(object){
	if(!assayDataValidMembers(assayData(object), c("alleleA", "alleleB", "call", "callProbability"))){
		return("assay data members must be 'alleleA', 'alleleB', 'call', 'callProbability'")
	}
	if(length(batch(object)) != ncol(object)){
 		return("'batch' must be the same length as the number of samples.  ")
 	}
	if(nrow(object) > 0){
		genome <- genomeBuild(object)
		if(!genome %in% c("hg18", "hg19")) print("Supported entries for genome are hg18 and hg19")
	}
	msg <- isValidGenomeAnnotatedDataFrame(featureData(object))
	if(is.character(msg)) return(msg) else TRUE
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


setMethod("initialize", "BeadStudioSet",
	  function(.Object,
		   assayData=assayDataNew(baf = baf, lrr = lrr, ...),
		   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
		   featureData = GenomeAnnotatedDataFrameFrom(assayData, annotation),
		   experimentData = new("MIAME"),
		   annotation = character(),
		   protocolData = phenoData[,integer(0)],
		   baf = new("matrix"),
		   lrr = matrix(numeric(),
                		   nrow=nrow(baf),
		                   ncol=ncol(baf),
                    		   dimnames=dimnames(baf)),
		   genome=c("hg19", "hg18"),
		   ...) {
		  .Object <- callNextMethod(.Object,
					    assayData = assayData,
					    phenoData = phenoData,
					    featureData = featureData,
					    experimentData = experimentData,
					    annotation = annotation,
					    protocolData = protocolData, ...)
		  if(nrow(assayData[["baf"]]) == 0){
			  .Object@genome <- ""
		  } else{
			  .Object@genome <- match.arg(genome)
		  }
	return(.Object)
})

setValidity("BeadStudioSet", function(object) {
	if(!is.null(lrr(object))){
		if(nrow(lrr(object)) > 0)
			if(!is.integer(lrr(object)[,1])) return("lrr should be a matrix of integers (original scale * 100). Use integerMatrix(x, 100) for converting 'x' to a matrix of integers.")
		if(nrow(baf(object)) > 0){
			b <- baf(object)[,1]
			if(!is.integer(b)) return("B allele frequencies should be a matrix of integers (original scale * 1000). See integerMatrix(x, 1000) for converting 'x' to a matrix of integers.")
		}
	}
	return(all(is.element(c("lrr","baf"), assayDataElementNames(object))))
})
