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
	if(!missing(end) && !missing(start))
		ranges <- IRanges(start, end)
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
	new("RangedDataCNV", ranges=ranges(rd), values=IRanges:::values(rd))
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
setMethod("coverage", signature(object="RangedDataCNV"), function(object) object$num.mark)
setMethod("sampleNames", signature(object="RangedDataCNV"), function(object) object$id)
setMethod("chromosome", signature(object="RangedDataCNV"), function(object) object$chrom)
setAs("RangedData", "RangedDataCBS", function(from){
	RangedDataCBS(ranges=ranges(from),
		      values=values(from))
})
setAs("RangedData", "RangedDataHMM", function(from){
	RangedDataHMM(ranges=ranges(from),
		      values=values(from))
})
