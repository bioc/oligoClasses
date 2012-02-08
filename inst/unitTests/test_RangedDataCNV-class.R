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
