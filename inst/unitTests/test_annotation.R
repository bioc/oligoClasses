test_annotation <- function(){
	library(oligoClasses)
	library(genomewidesnp6Crlmm)
	m <- matrix(NA, 2,1, dimnames=list(c("SNP_A-8575125", "CN_473963"), NULL))
	gad <- GenomeAnnotatedDataFrameFrom(m, "genomewidesnp6Crlmm")
	checkIdentical(position(gad), c(564621L, 61735L))
	checkIdentical(chromosome(gad), rep(1L,2))
}