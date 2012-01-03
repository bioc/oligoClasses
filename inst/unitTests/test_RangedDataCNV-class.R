make_test_RangedDataCNV <- function(){
	new("RangedDataCNV",
	    chrom=1:3,
	    coverage=rep(10L,3),
	    sampleId=letters[1:3])
}

test_RangedDataCNV_construction <- function(){
	checkException(RangedDataCNV(), silent=TRUE)
	checkTrue(validObject(RangedDataCNV()))
	checkTrue(validObject(new("RangedDataCNV")))
	checkIdentical(RangedDataCNV(chrom=1:3, coverage=rep(10L,3),
				     sampleId=letters[1:3]),
		       make_test_RangedDataCNV())
}
