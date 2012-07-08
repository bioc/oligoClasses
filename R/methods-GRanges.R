##GRangesHMM <- function(seqnames = Rle(), ranges = IRanges(),
##		       strand = Rle("*", length(seqnames)),
##		       state=integer(),
##		       numberProbes=integer(),
##		       ...,
##		       seqlengths = # this default is more accurate than in newGRanges
##		       structure(rep(NA_integer_, length(unique(seqnames))),
##				 names = levels(as.factor(runValue(as(seqnames, "Rle")))))){
##	if(is.numeric(state)){
##		state <- Rle(factor(state))
##	} else {
##		if(!is(state, "Rle")) stop("state must be Rle or numeric")
##	}
##	if (missing(seqlengths)) # avoid potentially expensive seqnames conversion
##		GenomicRanges:::newGRanges("GRangesHMM", seqnames = seqnames, ranges = ranges, strand = strand,
##					   state=state,
##					   numberProbes=as.integer(numberProbes),
##					   ...)
##	else GenomicRanges:::newGRanges("GRangesHMM", seqnames = seqnames, ranges = ranges,
##					strand = strand,
##					state=state,
##					numberProbes=as.integer(numberProbes),
##					..., seqlengths = seqlengths)
##}
##setMethod("initialize", "GRangesHMM",
##	  function(.Object, numberProbes=integer(), state=Rle(), elementMetadata=DataFrame(state=state, numberProbes=numberProbes), ...){
##		  callNextMethod(.Object, elementMetadata=elementMetadata, ...)
##	  })
##setGeneric("as.GRangesHMM", function(x, seqlengths, ...) standardGeneric("as.GRangesHMM"))
##setAs("GRanges", "GRangesHMM",
##      function(from, to){
##	      state <- elementMetadata(from)$state
##	      numberProbes <- elementMetadata(from)$numberProbes
##	      gr <- GRangesHMM(paste("chr", chromosome(from), sep=""),
##			       IRanges(start(from), end(from)),
##			       state=state(from),
##			       numberProbes=coverage2(from))
##      })
##GRangesHMMList <- function(...)
##{
##    listData <- list(...)
##    if (length(listData) == 0L) {
##        unlistData <- GRangesHMM()
##    } else {
##        if (length(listData) == 1L && is.list(listData[[1L]]))
##            listData <- listData[[1L]]
##        if (!all(sapply(listData, is, "GRangesHMM")))
##            stop("all elements in '...' must be GRangeHMM objects")
##        unlistData <- suppressWarnings(do.call("c", unname(listData)))
##    }
##    end <- cumsum(elementLengths(unname(listData)))
##    ans <- IRanges:::newCompressedList("GRangesHMMList",
##               unlistData,
##               end = end, NAMES = names(listData),
##               elementMetadata = new("DataFrame", nrows = length(listData)))
##    validObject(ans)
##    ans
##}

##---------------------------------------------------------------------------
##
## convenience functions for GRanges
##
##---------------------------------------------------------------------------
setMethod("chromosome", "GRanges", function(object) as.character(seqnames(object)))
setMethod("coverage2", "GRanges", function(object) elementMetadata(object)$numberProbes)
setMethod("numberProbes", "GRanges", function(object) elementMetadata(object)$numberProbes)
setMethod("sampleNames", "GRanges", function(object) as.character(elementMetadata(object)$sample))
setMethod("state", "GRanges", function(object) elementMetadata(object)$state)

##---------------------------------------------------------------------------
##
## convenience functions for GRangesList
##
##---------------------------------------------------------------------------
setMethod("chromosome", signature(object="GRangesList"), function(object) seqnames(object))
setMethod("coverage2", signature(object="GRangesList"), function(object) numberProbes(object))
setMethod("numberProbes", signature(object="GRangesList"), function(object){
	new2("CompressedRleList", unlistData=Rle(elementMetadata(object@unlistData)$numberProbes), partitioning=object@partitioning, check=FALSE)
})
setMethod("sampleNames", signature(object="GRangesList"), function(object) names(object))
setMethod("state", signature(object="GRangesList"), function(object){
	new2("CompressedRleList", unlistData=Rle(elementMetadata(object@unlistData)$state), partitioning=object@partitioning, check=FALSE)
})








