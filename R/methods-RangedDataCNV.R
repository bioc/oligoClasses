##setValidity("RangedDataCNV", function(object){
##	all(c("chrom", "id", "num.mark") %in% colnames(object))
##})
##setValidity("RangedDataCBS", function(object){
##	if(nrow(object) > 0){
##		all(c("seg.mean", "start.index", "end.index") %in% colnames(object))
##	}
##})
##setValidity("RangedDataHMM", function(object) "state" %in% colnames(object))
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
	.Defunct("RangedDataCNV is defunct. Use GenomicRanges instead.")
}

RangedDataCBS <- function(ranges=IRanges(),
			  seg.mean=vector("numeric", length(ranges)), ...){
	.Defunct("RangedDataCBS is defunct.  Use GRanges instead")
}

RangedDataHMM <- function(ranges=IRanges(),
			  state=vector("integer", length(ranges)), ...){
	.Defunct("RangedDataHMM is defunct. Use GRanges instead")
}

setMethod("state", signature(object="RangedDataCNV"), function(object) .Defunct())
setMethod("coverage2", signature(object="RangedDataCNV"), function(object) .Defunct())
setMethod("mean", signature(x="RangedDataCBS"), function(x,...) .Defunct())
setMethod("sampleNames", signature(object="RangedDataCNV"), function(object) .Defunct())
setMethod("chromosome", signature(object="RangedDataCNV"), function(object, na.rm=FALSE) .Defunct())
setMethod("findOverlaps", signature(query="RangedDataCNV", subject="SnpSet"),
	  function (query, subject, maxgap = 0L, minoverlap = 1L, type = c("any",
								  "start", "end", "within", "equal"), select = c("all", "first",
												      "last", "arbitrary"), ...){
		  findOverlaps(query=query, subject=featureData(subject),
			       maxgap=maxgap,
			       minoverlap=minoverlap,
			       type=type,
			       select=select, ...)
		  .Defunct("findOvelaps method for RangedDataCNV and SnpSet is defunct")
	  })

setMethod("findOverlaps", signature(query="RangedDataCNV", subject="CNSet"),
	  function (query, subject, maxgap = 0L, minoverlap = 1L, type = c("any",
								  "start", "end", "within", "equal"), select = c("all", "first",
												      "last", "arbitrary"), ...){
		  .Defunct("findOvelaps method for RangedDataCNV and CNSet is defunct")
	  })

setMethod("findOverlaps", signature(query="RangedDataCNV", subject="AnnotatedDataFrame"),
	  function (query, subject, maxgap = 0L, minoverlap = 1L, type = c("any",
								  "start", "end", "within", "equal"), select = c("all", "first",
												      "last", "arbitrary"), ...){
		  .Defunct("findOverlaps for RangedDataCNV and AnnotatedDataFrame is defunct")
	  })

setMethod("findOverlaps", signature(query="AnnotatedDataFrame", subject="RangedDataCNV"),
	  function (query, subject, maxgap = 0L, minoverlap = 1L, type = c("any",
								  "start", "end", "within", "equal"), select = c("all", "first",
												      "last", "arbitrary"), ...){
		  .Defunct("findOverlaps for AnnotatedDataFrame and RangedDataCNV is defunct")
	  })


setMethod("findOverlaps", signature(query="RangedDataCNV",
				    subject="RangedDataCNV"),
	  function(query, subject, maxgap = 0L, minoverlap = 1L,
		   type = c("any", "start", "end", "within", "equal"),
		   select = c("all", "first", "last", "arbitrary"), ...){
		  .Defunct("findOverlaps for RangedDataCNV is defunct")
	  })

setMethod("findOverlaps", signature(query="RangedDataHMM",
				    subject="RangedDataHMM"),
	  function(query, subject, maxgap = 0L, minoverlap = 1L,
		   type = c("any", "start", "end", "within", "equal"),
		   select = c("all", "first", "last", "arbitrary"), ...){
		  .Defunct("findOverlaps for RangedDataHMM is defunct")
	  })


setReplaceMethod("sampleNames", signature(object="RangedDataCNV",
					  value="character"),
		 function(object, value){
			 .Defunct("sampleNames<- defunct for RangedDataCNV")
		 })
setMethod("genomeBuild", signature(object="GRanges"),
	  function(object) metadata(object)[["genome"]])


##setAs("RangedDataHMM", "GRangesList", function(from, to){
##	GRangesListFromRangedDataHMM(from)
##})
##
##setAs("RangedDataCNV", "GRangesList", function(from, to){
##	GRangesListFromRangedDataCNV(from)
##})
##GRangesListFromRangedDataCNV <- function(object, build, ...){
##	index <- split(seq_len(nrow(object)), sampleNames(object))
##	notmissingbuild <- !missing(build)
##	if(notmissingbuild) sl <- getSequenceLengths(build)
##	grl <- vector("list", length(index))
##	for(i in seq_along(index)){
##		j <- index[[i]]
##		gr <- GRanges(paste("chr", chromosome(object)[j], sep=""),
##			      IRanges(start(object)[j], end(object)[j]))
##		if(notmissingbuild) seqlengths(gr) <- sl[match(unique(seqnames(gr)), names(sl))]
##		elementMetadata(gr)$numberProbes <- coverage2(object)[j]
##		grl[[i]] <- gr
##	}
##	grl <- GRangesList(grl)
##	names(grl) <- names(index)
##	grl
##}
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
