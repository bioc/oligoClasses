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
		  start <- start(query)
		  end <- end(query)
		  CHR <- chromosome(query)
		  ##featuresInXlim(object, start=start(range), end=end(range), CHR=range$chrom, ...)
##		  if("frame" %in% names(list(...))) {
##			  frame <- list(...)[["frame"]]
##		  } else frame <- rep(0, nrow(query))
##		  if(any(frame > 0)){
##			  data(chromosomeAnnotation, package="SNPchip")
##			  chr.end <- chromosomeAnnotation[CHR, "chromosomeSize"]
##			  start <- start-frame
##			  start[start < 0] <- 0
##			  end <- end+frame
##			  end[end > chr.end] <- chr.end
##		  }
		  ir.query <- IRanges(start, end)
		  ## depends on platform
		  ir.subject <- IRanges(position(subject)-12, position(subject)+12)
		  res <- findOverlaps(query=ir.query,
				      subject=ir.subject,
				      maxgap=maxgap,
				      minoverlap=minoverlap,
				      type=type,
				      select=select, ...)
		  mm <- matchMatrix(res)
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
		  }
		  res@matchMatrix <- mm
		  return(res)
	  })


setMethod("findOverlaps", signature(query="RangedDataCNV",
				    subject="RangedDataCNV"),
	  function(query, subject, maxgap = 0L, minoverlap = 1L,
		   type = c("any", "start", "end", "within", "equal"),
		   select = c("all", "first", "last", "arbitrary"),
		   match.id=TRUE, ...){
		  irq <- IRanges(start(query), end(query))
		  irs <- IRanges(start(subject), end(subject))
		  chrq <- chromosome(query)
		  chrs <- chromosome(subject)
		  if("match.id" %in% names(list(...))){
			  match.id <- list(...)[["match.id"]]
		  } else match.id <- TRUE
		  if(match.id){
			  idq <- sampleNames(query)
			  ids <- sampleNames(subject)
		  }
		  res <- findOverlaps(irq, irs, maxgap=maxgap,
				      minoverlap=minoverlap,
				      type=type,
				      select=select,...)
		  mm <- matchMatrix(res)
		  chrq <- chrq[mm[,1]]
		  chrs <- chrs[mm[,2]]
		  if(match.id){
			  idq <- idq[mm[,1]]
			  ids <- ids[mm[,2]]
			  index <- chrq == chrs & idq == ids
		  } else index <- chrq == chrs
		  res@matchMatrix <- mm[index, , drop=FALSE]
		  return(res)
	  })

setMethod("findOverlaps", signature(query="RangedDataHMM",
				    subject="RangedDataHMM"),
	  function(query, subject, maxgap = 0L, minoverlap = 1L,
		   type = c("any", "start", "end", "within", "equal"),
		   select = c("all", "first", "last", "arbitrary"), ...){
		  stateq <- state(query)
		  states <- state(subject)
		  res <- callNextMethod(...)
		  mm <- matchMatrix(res)
		  stateq <- stateq[mm[, 1]]
		  states <- states[mm[, 2]]
		  index <- which(stateq == states)
		  res@matchMatrix <- mm[index, , drop=FALSE]
		  return(res)
	  })

setReplaceMethod("sampleNames", signature(object="RangedDataCNV",
					  value="character"),
		 function(object, value){
			 object$id <- value
			 return(object)
		 })
