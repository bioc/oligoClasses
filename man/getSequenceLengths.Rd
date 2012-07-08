\name{getSequenceLengths}
\alias{getSequenceLengths}
%\alias{setSequenceLengths}

\title{Load chromosome sequence lengths for UCSC genome build hg18 or hg19}
\description{Load chromosome sequence lengths for UCSC genome build hg18 or hg19}
\usage{
getSequenceLengths(build)
%setSequenceLengths(build, names)
}
\arguments{
  \item{build}{
    character string: "hg18" or "hg19"
  }

%  \item{names}{
%
%    Chromosome names  (e.g, 'chr1', 'chr2', ...)
%
%  }

}
\details{The chromosome sequence lengths for UCSC builds hg18 and hg19
  were extracted from the packages BSgenome.Hsapiens.UCSC.hg18 and
  BSgenome.Hsapiens.UCSC.hg19, respectively.
}
\value{

  Names integer vector of chromosome lengths.

}

\author{
R. Scharpf
}

\examples{
getSequenceLengths("hg18")
getSequenceLengths("hg19")

if(require("GenomicRanges")){
## from GenomicRanges
	sl <- getSequenceLengths("hg18")[c("chr1", "chr2", "chr3")]
	gr <-
		GRanges(seqnames =
			Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
			ranges =
			IRanges(1:10, width = 10:1, names = head(letters,10)),
			strand =
			Rle(strand(c("-", "+", "*", "+", "-")),
			    c(1, 2, 2, 3, 2)),
			score = 1:10,
			GC = seq(1, 0, length=10),
			seqlengths=sl)
	metadata(gr) <- list(genome="hg18")
	gr
	metadata(gr)
}


}
\keyword{attribute}

