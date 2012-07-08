##make_test_RangedDataCNV <- function(){
##	new("RangedDataCNV",
##	    start=1:3, end=2:4,
##	    chromosome=1:3,
##	    coverage=rep(10L,3),
##	    sampleId=letters[1:3])
##}

test_RangedDataCNV_construction <- function(){
	##checkException(RangedDataCNV(), silent=TRUE)
	checkTrue(validObject(RangedDataCNV(start=1:3, end=2:4,
					    chrom=1:3,
					    coverage=rep(10L,3),
					    sampleId=letters[1:3])))
}

test_findOverlaps <- function(){
	if(require(IRanges)){
		data(oligoSetExample)
		## RangedDataHMM and AnnotatedDataFrame/GenomeAnnotatedDataFrame
		rd.hmm <- RangedDataHMM(IRanges(start=49825223, end=54637167),
					chrom=1, state=3)
		index <- subjectHits(findOverlaps(rd.hmm, featureData(oligoSet)))
		checkIdentical(index, seq(997, 1098))

		## swap query and subject
		index2 <- queryHits(findOverlaps(featureData(oligoSet), rd.hmm))
		checkIdentical(index, index2)
	}
}

test_coersionToGRanges <- function(){
	library(GenomicRanges)
	data(oligoSetExample)
	## RangedDataHMM and AnnotatedDataFrame/GenomeAnnotatedDataFrame
	rd.hmm <- RangedDataHMM(IRanges(start=49825223, end=54637167),
				chrom=1, state=3)
	if(FALSE){ ## not exported yet
		checkTrue(validObject(new("GRangesHMM")))
		gr <- GRangesHMM(seqnames="chr1", ranges=IRanges(1,3),
				 strand="*", seqlengths=c("chr1"=1000),
				 numberProbes=5L, state=3L)
		checkTrue(validObject(gr))
		library(BSgenome.Hsapiens.UCSC.hg18)
		library(BSgenome.Hsapiens.UCSC.hg19)
		seqlengths <- seqlengths(Hsapiens)
		obj <- as.GRangesHMM(rd.hmm, seqlengths=seqlengths(Hsapiens))
		checkTrue(validObject(obj))
		checkException(as.GRangesHMM(rd.hmm))
	}
}
