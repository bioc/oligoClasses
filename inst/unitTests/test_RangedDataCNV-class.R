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
