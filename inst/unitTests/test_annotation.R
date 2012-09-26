test_annotation <- function(){
	library(oligoClasses)
	library(genomewidesnp6Crlmm)
	m <- matrix(NA, 2,1, dimnames=list(c("SNP_A-8575125", "CN_473963"), NULL))
	gad <- GenomeAnnotatedDataFrameFrom(m, "genomewidesnp6Crlmm", genome="hg19")
	##checkIdentical(position(gad), c(564621L, 61723L))
	checkIdentical(chromosome(gad), rep(1L,2))

	## genome annotated data frame with additional columns
	gd <- new("GenomeAnnotatedDataFrame",
		  position=1:10,
		  chrom=1:10,
		  isSnp=TRUE,
		  arm=1)
	checkTrue(validObject(gd))
	## other approach
	library(Biobase)
	pD <- pData(gd)
	mD <- varMetadata(gd)
	gd <- new("GenomeAnnotatedDataFrame",
		  data=pD,
		  varMetadata=mD)
	checkTrue(validObject(gd))

	gd <- gd[1:5, ]
	checkTrue(validObject(gd))

	x <- matrix(1:25, 5, 5, dimnames=list(c("rs10000092","rs1000055", "rs100016", "rs10003241", "rs10004197"), NULL))
	## preferred
	gd <- GenomeAnnotatedDataFrameFrom(x, annotationPkg="human370v1cCrlmm", genome="hg18")
	checkTrue(is(gd, "GenomeAnnotatedDataFrame"))
	## searches for hg19 by default
	## checkException(gd <- GenomeAnnotatedDataFrameFrom(x, annotationPkg="human370v1cCrlmm"))
	checkTrue(is(gd, "GenomeAnnotatedDataFrame"))
	gd <- GenomeAnnotatedDataFrameFrom(x, annotationPkg="human370v1cCrlmm", genome="")
	checkTrue(is(gd, "GenomeAnnotatedDataFrame"))
	## request a build that is not available
	## checkException(GenomeAnnotatedDataFrameFrom(x, annotationPkg="human370v1cCrlmm", genome="hg19"))
}


