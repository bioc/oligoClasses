test_conversions <- function(){
	p <- matrix(runif(20), nc=2)
	integerRepresentation <- as.integer(-1000*log(1-p))
	int2 <- p2i(p)
	checkTrue(all.equal(integerRepresentation, int2))
}

test_oligoSnpSet <- function(){
	data(oligoSetExample)
	checkTrue(validObject(as(oligoSet, "SnpSet2")))
}

test_makeFeatureRanges <- function(){
	data(oligoSetExample)
	gr <- makeFeatureGRanges(featureData(oligoSet), genome=genomeBuild(oligoSet))
	checkTrue(validObject(gr))
	gr2 <- makeFeatureGRanges(oligoSet)
	checkIdentical(gr, gr2)
}

##test_RangedDataHMM2GRanges <- function(){
##	if(require(VanillaICE)){
##		data(hmmResults, package="VanillaICE")
##		checkTrue(validObject(as(hmmResults, "GRanges")))
##		obj <- as(hmmResults, "GRangesList")
##		checkTrue(validObject(obj))
##		checkEquals(names(obj), unique(sampleNames(hmmResults)))
##	}
##}
