setValidity("RangedDataCNV", function(object){
	all(c("chrom", "id", "num.mark") %in% colnames(object))
})

setValidity("RangedDataCBS", function(object){
	if(nrow(object) > 0){
		all(c("seg.mean", "start.index", "end.index") %in% colnames(object))
	}
})
setValidity("RangedDataHMM", function(object) "state" %in% colnames(object))

RangedDataCNV <- function(ranges=IRanges(),
			  values,
			  start,
			  end,
			  chromosome,
			  coverage,
			  sampleId,
			  startIndexInChromosome,
			  endIndexInChromosome,
			  ...){
	if(!missing(ranges) & !missing(values)){
		object <- new("RangedDataCNV", ranges=ranges, values=values)
		return(object)
	}
	if(!missing(end) && !missing(start)){
		ranges <- IRanges(start, end)
	}
	if(missing(chromosome))
		chromosome <- vector("integer", length(ranges))
	if(missing(coverage))
		coverage <- vector("integer", length(ranges))
	if(missing(sampleId))
		sampleId <- vector("character", length(ranges))
	if(missing(startIndexInChromosome))
		startIndexInChromosome <- vector("integer", length(ranges))
	if(missing(endIndexInChromosome))
		endIndexInChromosome <- vector("integer", length(ranges))
	rd <- RangedData(ranges,
			 chrom=chromosome,
			 num.mark=coverage,
			 id=sampleId,
			 start.index=startIndexInChromosome,
			 end.index=endIndexInChromosome, ...)##, ...)
	if(nrow(rd) > 0){
		obj <- new("RangedDataCNV", ranges=ranges(rd), values=IRanges:::values(rd))
	} else {
		obj <- new("RangedDataCNV")
	}
	return(obj)
}

RangedDataCBS <- function(ranges=IRanges(),
			  seg.mean=vector("numeric", length(ranges)), ...){
	rd <- RangedDataCNV(ranges=ranges, seg.mean=seg.mean, ...)
	new("RangedDataCBS", ranges=ranges(rd), values=values(rd))
}

RangedDataHMM <- function(ranges=IRanges(),
			  state=vector("integer", length(ranges)), ...){
	rd <- RangedDataCNV(ranges=ranges, state=state, ...)
	new("RangedDataHMM", ranges=ranges(rd), values=values(rd))
}

setMethod("state", signature(object="RangedDataCNV"), function(object) object$state)
##setMethod("nMarkers", signature(object="RangedDataCNV"), function(object) object$num.mark)
##setMethod("coverage", signature(object="RangedDataCNV"), function(object) coverage2(object))
setMethod("coverage2", signature(object="RangedDataCNV"), function(object) object$num.mark)
setMethod("mean", signature(x="RangedDataCBS"), function(x,...) x$seg.mean)
setMethod("sampleNames", signature(object="RangedDataCNV"), function(object) object$id)
setMethod("chromosome", signature(object="RangedDataCNV"), function(object, na.rm=FALSE) object$chrom)
setAs("RangedData", "RangedDataCBS", function(from){
	RangedDataCBS(ranges=ranges(from),
		      values=values(from))
})
setAs("RangedData", "RangedDataHMM", function(from){
	RangedDataHMM(ranges=ranges(from),
		      values=values(from))
})

setMethod("findOverlaps", signature(query="RangedDataCNV", subject="SnpSet"),
	  function (query, subject, maxgap = 0L, minoverlap = 1L, type = c("any",
								  "start", "end", "within", "equal"), select = c("all", "first",
												      "last", "arbitrary"), ...){
		  findOverlaps(query=query, subject=featureData(subject),
			       maxgap=maxgap,
			       minoverlap=minoverlap,
			       type=type,
			       select=select, ...)
	  })

setMethod("findOverlaps", signature(query="RangedDataCNV", subject="CNSet"),
	  function (query, subject, maxgap = 0L, minoverlap = 1L, type = c("any",
								  "start", "end", "within", "equal"), select = c("all", "first",
												      "last", "arbitrary"), ...){
		  findOverlaps(query=query, subject=featureData(subject),
			       maxgap=maxgap,
			       minoverlap=minoverlap,
			       type=type,
			       select=select, ...)
	  })

setMethod("findOverlaps", signature(query="RangedDataCNV", subject="AnnotatedDataFrame"),
	  function (query, subject, maxgap = 0L, minoverlap = 1L, type = c("any",
								  "start", "end", "within", "equal"), select = c("all", "first",
												      "last", "arbitrary"), ...){
		  subject2 <- subject
		  nachrom <- is.na(chromosome(subject)) | chromosome(subject) > 24
		  if(any(nachrom)){
			  subject <- subject[!nachrom, ]
		  }
		  CHR <- chromosome(query)
		  ir.query <- IRanges(start(query), end(query))
		  ir.subject <- IRanges(position(subject), position(subject))
		  res <- findOverlaps(query=ir.query,
				      subject=ir.subject,
				      maxgap=maxgap,
				      minoverlap=minoverlap,
				      type=type,
				      select=select, ...)
		  mm <- as.matrix(res)
		  ## remove matches that are not the same chromosome
		  subj.index <- mm[,2]
		  quer.index <- mm[, 1]
		  same.chrom <- chromosome(query)[quer.index] == chromosome(subject)[subj.index]
		  ##subj.index <- subj.index[same.chrom]
		  ##quer.index <- quer.index[same.chrom]
		  ## narrow the hits to those that are in the same chromosome
		  ##which(position(object) >= start & position(object) <= end & chromosome(object) == CHR)
		  ##} else which(pos >= start & pos <= end & chrom == CHR)
		  mm <- mm[same.chrom, , drop=FALSE]
		  res <- res[same.chrom, ]
		  ## Now, map the subject indices back to the indices in
		  ## the original object
		  if(any(nachrom)){
			  ## 'sampleNames' can be used to retrieve the feature names
			  subject.index <- mm[, 2]
			  fns.subject <- sampleNames(subject)[subject.index]
			  fns.subject2 <- sampleNames(subject2)
			  index <- match(fns.subject, fns.subject2)
			  stopifnot(all(!is.na(index)))
			  mm[, 2] <- index
			  res <- res[index, ]
		  }
		  ##res@as.matrix <- mm
		  return(res)
	  })

setMethod("findOverlaps", signature(query="AnnotatedDataFrame", subject="RangedDataCNV"),
	  function (query, subject, maxgap = 0L, minoverlap = 1L, type = c("any",
								  "start", "end", "within", "equal"), select = c("all", "first",
												      "last", "arbitrary"), ...){
		  query2 <- query
		  nachrom <- is.na(chromosome(query)) | chromosome(query) > 24
		  if(any(nachrom)){
			  query <- query[!nachrom, ]
		  }
		  ##start <- start(subject)
		  ##end <- end(subject)
		  CHR <- chromosome(subject)
		  ir.subject <- IRanges(start(subject), end(subject))
		  ir.query <- IRanges(position(query), position(query))
		  res <- findOverlaps(query=ir.query,
				      subject=ir.subject,
				      maxgap=maxgap,
				      minoverlap=minoverlap,
				      type=type,
				      select=select, ...)
		  mm <- as.matrix(res)
		  subj.index <- mm[,2]
		  quer.index <- mm[, 1]
		  same.chrom <- chromosome(query)[quer.index] == chromosome(subject)[subj.index]
		  mm <- mm[same.chrom, , drop=FALSE]
		  res <- res[same.chrom, ]
		  if(any(nachrom)){
			  query.index <- mm[, 1]
			  fns.query <- sampleNames(query)[query.index]
			  fns.query2 <- sampleNames(query2)
			  index <- match(fns.query, fns.query2)
			  stopifnot(all(!is.na(index)))
			  mm[, 1] <- index
			  res <- res[index, ]
		  }
		  ##res@as.matrix <- mm
		  return(res)
	  })

findOverlapsForRangedDataCNV <- function(query, subject,
					 maxgap = 0L, minoverlap = 1L,
					 type = c("any", "start", "end", "within", "equal"),
					 select = c("all", "first", "last", "arbitrary"),
					 ...){
	##ids are irrelevant
	irq <- IRanges(start(query), end(query))
	irs <- IRanges(start(subject), end(subject))
	chrq <- chromosome(query)
	chrs <- chromosome(subject)
	res <- findOverlaps(irq, irs, maxgap=maxgap,
			    minoverlap=minoverlap,
			    type=type,
			    select=select, ...)
	mm <- as.matrix(res)
	chrq <- chrq[mm[,1]]
	chrs <- chrs[mm[,2]]
	idsq <- sampleNames(query)[mm[, 1]]
	idss <- sampleNames(subject)[mm[,2]]
	index <- which(chrq == chrs & idsq == idss)
	res <- res[index, ]
	return(res)
}


setMethod("findOverlaps", signature(query="RangedDataCNV",
				    subject="RangedDataCNV"),
	  function(query, subject, maxgap = 0L, minoverlap = 1L,
		   type = c("any", "start", "end", "within", "equal"),
		   select = c("all", "first", "last", "arbitrary"), ...){
		  findOverlapsForRangedDataCNV(query=query, subject=subject,
					       maxgap=maxgap,
					       type=type,
					       select=select,
					       ...)
	  })

##findOverlapsForRangedDataWithIds <- function(query, subject,
##					     maxgap = 0L, minoverlap = 1L,
##					     type = c("any", "start", "end", "within", "equal"),
##					     select = c("all", "first", "last", "arbitrary"),
##					     ...){
##	res <- findOverlapsForRangedDataCNV(query, subject,
##					    maxgap = maxgap,
##					    minoverlap = minoverlap,
##					    type = type,
##					    select = select,
##					    ...)
##	mm <- as.matrix(res)
##	idsq <- sampleNames(query)[mm[, 1]]
##	idss <- sampleNames(subject)[mm[,2]]
##	index <- which(idsq==idss)
##	res[index, ]
##}

setMethod("findOverlaps", signature(query="RangedDataHMM",
				    subject="RangedDataHMM"),
	  function(query, subject, maxgap = 0L, minoverlap = 1L,
		   type = c("any", "start", "end", "within", "equal"),
		   select = c("all", "first", "last", "arbitrary"), ...){
		  res <- findOverlapsForRangedDataCNV(query,
						      subject,
						      maxgap=maxgap,
						      minoverlap=minoverlap,
						      type=type,
						      select=select,
						      ...)
		  mm <- as.matrix(res)
		  stateq <- state(query)
		  states <- state(subject)
		  ##res <- callNextMethod(...)
		  stateq <- stateq[mm[, 1]]
		  states <- states[mm[, 2]]
		  index <- which(stateq == states)
		  res <- res[index, ]
		  return(res)
	  })


setReplaceMethod("sampleNames", signature(object="RangedDataCNV",
					  value="character"),
		 function(object, value){
			 object$id <- value
			 return(object)
		 })



##stack2 <- function(x){
##	ids <- rep(names(x), elementLengths(x))
##	x <- stack(x)
##	elementMetadata(x)$id <- ids
##	x
##}

setMethod("makeFeatureGRanges", signature(object="gSet"),
	  function(object, ...){
		  sl <- getSequenceLengths(genomeBuild(object))
		  gr <- GRanges(paste("chr", chromosome(object), sep=""),
				IRanges(position(object), width=1))
		  seqlengths(gr) <- sl[match(unique(seqnames(gr)), names(sl))]
		  return(gr)
	  })

setMethod("genomeBuild", signature(object="GRanges"),
	  function(object) metadata(object)[["genome"]])

setMethod("makeFeatureGRanges", signature(object="GenomeAnnotatedDataFrame"),
	  function(object, genome, ...){
		  sl <- getSequenceLengths(genome)
		  chrom <- integer2chromosome(chromosome(object))
		  gr <- GRanges(paste("chr", chrom, sep=""),
				IRanges(position(object), width=1))
		  nms <- names(seqlengths(gr))
		  sl <- sl[match(nms, names(sl))]
		  seqlengths(gr) <- sl
		  return(gr)
	  })


GRangesListFromRangedDataHMM <- function(object, build, ...){
	index <- split(seq_len(nrow(object)), sampleNames(object))
	notmissingbuild <- !missing(build)
	if(notmissingbuild) sl <- getSequenceLengths(build)
	grl <- vector("list", length(index))
	for(i in seq_along(index)){
		j <- index[[i]]
		gr <- GRanges(paste("chr", chromosome(object)[j], sep=""),
			      IRanges(start(object)[j], end(object)[j]))
		if(notmissingbuild) seqlengths(gr) <- sl[match(unique(seqnames(gr)), names(sl))]
		elementMetadata(gr)$state <- state(object)[j]
		elementMetadata(gr)$numberProbes <- coverage2(object)[j]
		elementMetadata(gr)$LLR <- object$LLR[j]
		grl[[i]] <- gr
	}
	grl <- GRangesList(grl)
	names(grl) <- names(index)
	grl
}

setAs("RangedDataHMM", "GRangesList", function(from, to){
	GRangesListFromRangedDataHMM(from)
})

setAs("RangedDataCNV", "GRangesList", function(from, to){
	GRangesListFromRangedDataCNV(from)
})



GRangesListFromRangedDataCNV <- function(object, build, ...){
	index <- split(seq_len(nrow(object)), sampleNames(object))
	notmissingbuild <- !missing(build)
	if(notmissingbuild) sl <- getSequenceLengths(build)
	grl <- vector("list", length(index))
	for(i in seq_along(index)){
		j <- index[[i]]
		gr <- GRanges(paste("chr", chromosome(object)[j], sep=""),
			      IRanges(start(object)[j], end(object)[j]))
		if(notmissingbuild) seqlengths(gr) <- sl[match(unique(seqnames(gr)), names(sl))]
		elementMetadata(gr)$numberProbes <- coverage2(object)[j]
		grl[[i]] <- gr
	}
	grl <- GRangesList(grl)
	names(grl) <- names(index)
	grl
}

setMethod("findOverlaps", signature(query="GRangesList", subject="gSet"),
	  function (query, subject,
		    maxgap = 0L, minoverlap = 1L,
		    type = c("any", "start", "end", "within", "equal"),
		    select = c("all", "first", "last", "arbitrary"), ...){
		  frange <- makeFeatureGRanges(subject)
		  findOverlaps(query, frange, maxgap=maxgap,
			       minoverlap=minoverlap,
			       type=match.arg(type),
			       select=match.arg(select),
			       ...)
	  })
setMethod("findOverlaps", signature(query="GRanges", subject="gSet"),
	  function (query, subject,
		    maxgap = 0L, minoverlap = 1L,
		    type = c("any", "start", "end", "within", "equal"),
		    select = c("all", "first", "last", "arbitrary"), ...){
		  frange <- makeFeatureGRanges(subject)
		  findOverlaps(query, frange, maxgap=maxgap,
			       minoverlap=minoverlap,
			       type=match.arg(type),
			       select=match.arg(select),
			       ...)
	  })

coerceToGRanges <- function(range, build="hg18"){
	##chrlevels <- paste("chr", 1:22, sep="")
	chrlevels <- names(getSequenceLengths(build))
	chrom <- paste("chr", chromosome(range), sep="")
	chrlevels <- chrlevels[chrlevels %in% chrom]
	if(is(range, "RangedDataHMM")){
		gr <- GRanges(factor(chrom, levels=chrlevels),
			      IRanges(start(range), end(range)),
			      sample=sampleNames(range),
			      state=range$state,
			      numberProbes=coverage2(range),
			      seqlengths=setSequenceLengths(build, names=chrlevels))
	}
	if(is(range, "RangedDataCBS")){
		gr <- GRanges(factor(chrom, levels=chrlevels),
			      IRanges(start(range), end(range)),
			      sample=sampleNames(range),
			      seg.mean=range$seg.mean,
			      numberProbes=coverage2(range),
			      seqlengths=setSequenceLengths(build, names=chrlevels))
	}
	metadata(gr) <- list(genome=build)
	gr
	##sort(gr)
}
