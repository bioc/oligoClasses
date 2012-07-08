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
