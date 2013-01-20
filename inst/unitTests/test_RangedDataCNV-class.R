##test_coersionToGRanges <- function(){
##	library(GenomicRanges)
##	data(oligoSetExample)
##	## RangedDataHMM and AnnotatedDataFrame/GenomeAnnotatedDataFrame
##	rd.hmm <- RangedDataHMM(IRanges(start=49825223, end=54637167),
##				chrom=1, state=3)
##	if(FALSE){ ## not exported yet
##		checkTrue(validObject(new("GRangesHMM")))
##		gr <- GRangesHMM(seqnames="chr1", ranges=IRanges(1,3),
##				 strand="*", seqlengths=c("chr1"=1000),
##				 numberProbes=5L, state=3L)
##		checkTrue(validObject(gr))
##		library(BSgenome.Hsapiens.UCSC.hg18)
##		library(BSgenome.Hsapiens.UCSC.hg19)
##		seqlengths <- seqlengths(Hsapiens)
##		obj <- as.GRangesHMM(rd.hmm, seqlengths=seqlengths(Hsapiens))
##		checkTrue(validObject(obj))
##		checkException(as.GRangesHMM(rd.hmm))
##	}
##}

