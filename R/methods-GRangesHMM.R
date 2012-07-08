GRangesHMM <- function(seqnames = Rle(), ranges = IRanges(),
		       strand = Rle("*", length(seqnames)),
		       state=integer(),
		       numberProbes=integer(),
		       ...,
		       seqlengths = # this default is more accurate than in newGRanges
		       structure(rep(NA_integer_, length(unique(seqnames))),
				 names = levels(as.factor(runValue(as(seqnames, "Rle")))))){
	if(is.numeric(state)){
		state <- Rle(factor(state))
	} else {
		if(!is(state, "Rle")) stop("state must be Rle or numeric")
	}
	if (missing(seqlengths)) # avoid potentially expensive seqnames conversion
		GenomicRanges:::newGRanges("GRangesHMM", seqnames = seqnames, ranges = ranges, strand = strand,
					   state=state,
					   numberProbes=as.integer(numberProbes),
					   ...)
	else GenomicRanges:::newGRanges("GRangesHMM", seqnames = seqnames, ranges = ranges,
					strand = strand,
					state=state,
					numberProbes=as.integer(numberProbes),
					..., seqlengths = seqlengths)
}
setMethod("initialize", "GRangesHMM",
	  function(.Object, numberProbes=integer(), state=Rle(), elementMetadata=DataFrame(state=state, numberProbes=numberProbes), ...){
		  callNextMethod(.Object, elementMetadata=elementMetadata, ...)
	  })
setMethod("chromosome", "GRanges", function(object) seqnames(object))
setMethod("state", "GRanges", function(object) elementMetadata(object)$state)
setMethod("coverage2", "GRanges", function(object) elementMetadata(object)$numberProbes)
setGeneric("numberProbes", function(object) standardGeneric("numberProbes"))
setMethod("numberProbes", "GRanges", function(object) elementMetadata(object)$numberProbes)

setGeneric("as.GRangesHMM", function(x, seqlengths, ...) standardGeneric("as.GRangesHMM"))
setMethod("as.GRangesHMM", "RangedDataHMM",
	  function(x, seqlengths, ...)
  {
	  if (missing(seqlengths))
		  stop("must supply seqlengths")
	  GRangesHMM(seqnames=Rle(factor(paste("chr", chromosome(x), sep=""))),
		     IRanges(start(x), end(x)),
		     state=state(x),
		     numberProbes=coverage2(x),
		     seqlengths=seqlengths)
  })

##setAs("GRanges", "GRangesHMM",
##      function(from, to){
##	      state <- elementMetadata(from)$state
##	      numberProbes <- elementMetadata(from)$numberProbes
##	      gr <- GRangesHMM(paste("chr", chromosome(from), sep=""),
##			       IRanges(start(from), end(from)),
##			       state=state(from),
##			       numberProbes=coverage2(from))
##      })



