test_dataExamples <- function(){
	data(oligoSetExample)
	checkTrue(validObject(oligoSet))
}

test_GenomeAnnotatedDataFrame_construction <- function(){
	checkTrue(validObject(new("GenomeAnnotatedDataFrame")))
	checkTrue(validObject(GenomeAnnotatedDataFrameFrom(NULL)))
	data(locusLevelData)
	require(pd.mapping50k.hind240)


	tmp <- GenomeAnnotatedDataFrameFrom(locusLevelData[["genotypes"]],
					    annotationPkg=locusLevelData[["platform"]])
	checkTrue(validObject(tmp))
}

test_oligoSnpSet_construction <- function(){
	checkTrue(validObject(new("oligoSnpSet")))
	data(locusLevelData)
	require(pd.mapping50k.hind240)
	require(pd.mapping50k.xba240)
	oligoSet <- new("oligoSnpSet",
			copyNumber=log2(locusLevelData[["copynumber"]]/100),
			call=locusLevelData[["genotypes"]],
			callProbability=locusLevelData[["crlmmConfidence"]],
			annotation=locusLevelData[["platform"]])
	checkTrue(validObject(oligoSet))
	## instantiate oligoSnpSet with 0-row featureData
	oligoSet <- new("oligoSnpSet",
			copyNumber=log2(locusLevelData[["copynumber"]]/100),
			call=locusLevelData[["genotypes"]],
			callProbability=locusLevelData[["crlmmConfidence"]])
	checkTrue(validObject(oligoSet))
}

test_CopyNumberSet_construction <- function(){
	checkTrue(validObject(new("CopyNumberSet")))
	require(pd.mapping50k.hind240)
	require(pd.mapping50k.xba240)
	data(locusLevelData)
	cnset <- new("CopyNumberSet",
		     copyNumber=log2(locusLevelData[["copynumber"]]/100),
		     annotation=locusLevelData[["platform"]])
	checkTrue(validObject(cnset))
	## instantiate oligoSnpSet with 0-row featureData
	cnset <- new("CopyNumberSet",
		     copyNumber=log2(locusLevelData[["copynumber"]]/100))
	checkTrue(validObject(cnset))
}

test_CNSet_construction <- function(){
	checkTrue(validObject(new("CNSet")))
	a <- matrix(1:25, 5, 5, dimnames=list(letters[1:5], LETTERS[1:5]))
	tmp <- new("CNSet", alleleA=a, batch=rep("a", 5))
	checkTrue(validObject(tmp))
	tmp2 <- tmp[1:3, 1:2]
	checkTrue(validObject(tmp2))

	tmp2 <- new("CNSet", alleleA=a, batch=c(rep("a", 3), "b", "b"))
	checkTrue(validObject(tmp2))
	checkTrue(identical(batchNames(tmp2), c("a", "b", "grandMean")))


	require("genomewidesnp6Crlmm")
	fns <- c("SNP_A-2131660", "SNP_A-1967418", "SNP_A-1969580", "SNP_A-4263484",
		 "SNP_A-1978185", "SNP_A-4264431", "SNP_A-1980898", "SNP_A-1983139",
		 "SNP_A-4265735", "SNP_A-1995832")
	theCalls <- matrix(2, nc=2, nrow=10)
	A <- matrix(sample(1:1000, 20), 10,2)
	B <- matrix(sample(1:1000, 20), 10,2)
	p <- matrix(runif(20), nc=2)
	theConfs <- round(-1000*log2(1-p))
	rownames(A) <- rownames(B) <- rownames(theConfs) <- fns
	batch <- rep("a", ncol(A))
	obj <- new("CNSet",
		   alleleA=A,
		   alleleB=B,
		   call=theCalls,
		   callProbability=theConfs,
		   batch=batch,
		   annotation="genomewidesnp6")
	checkTrue(validObject(obj))
	checkTrue(identical(sampleNames(batchStatistics(obj)), batchNames(obj)))
	checkTrue(!is.null(batchNames(obj)))
	checkTrue(all(chromosome(obj) == 1))
}

test_GenomeAnnotatedDataFrameWithFF <- function(){
	## test instantiation from an object of class ff_matrix
	data(oligoSetExample)
	fdFromMatrix <- GenomeAnnotatedDataFrameFrom(locusLevelData[["genotypes"]],
						     annotationPkg=locusLevelData[["platform"]])

	if(require(ff)){
		ldPath(tempdir())
		gtMatrix <- locusLevelData[["genotypes"]]
		gts <- initializeBigMatrix(name="genotypes", initdata=gtMatrix, nr=nrow(gtMatrix), nc=ncol(gtMatrix), vmode="integer")
		rownames(gts) <- rownames(gtMatrix)
		fdFromFF <- GenomeAnnotatedDataFrameFrom(gts,
							 annotationPkg=locusLevelData[["platform"]])
		checkTrue(identical(fdFromMatrix, fdFromFF))
	}
}

test_BeadStudioSet <- function(){
	checkTrue(validObject(new("BeadStudioSet")))
}
