\name{GRanges-methods}
\docType{methods}
%\alias{coverage2,GRanges-method}
%\alias{coverage2,GRangesList-method}
\alias{state}
\alias{coverage2}
\alias{numberProbes}
\alias{coverage2,GRanges-method}
\alias{coverage2,GRangesList-method}
\alias{numberProbes,GRanges-method}
\alias{numberProbes,GRangesList-method}
\alias{findOverlaps,GRanges,gSet-method}
\alias{findOverlaps,GRangesList,gSet-method}
\alias{genomeBuild,GRanges-method}
\alias{sampleNames,GRanges-method}
\alias{sampleNames,GRangesList-method}
\alias{state,GRanges-method}
\alias{state,GRangesList-method}

\title{Methods for GRanges objects}

\description{Methods for GRanges objects}

\section{findOverlaps methods}{

  \describe{

    \item{}{
      \code{findOverlaps(query, subject, ...)}:

      Find the feature indices in \code{subject} that overlap the
      genomic intervals in \code{query}, where \code{query} is a
      \code{GRanges} object and subject is a \code{gSet}-derived object.
      Additional arguments to the \code{findOverlaps} method in the
      package \pkg{IRanges} can be passed through the \code{\dots}
      operator.

  }

}
}

\section{Accessors}{

  \code{object}  is an instance of the \code{GRanges} class.

  \describe{

    \item{}{
      \code{coverage2(object)}:

      For the \code{GRanges} and \code{GRangesList} objects returned by
      the hidden Markov model implemented in the "VanillaICE" package
      and the segmentation algorithm in the "MinimumDistance" package,
      the intervals are annotated by the number of probes (markers) for
      SNPs and nonpolymorphic regions.  \code{coverage2} and
      \code{numberProbes} are convenient accessors for these
      annotations.

    }

    \item{}{
      \code{genomeBuild(object)}:

      Accessor for the UCSC genome build.

    }

    \item{}{
      \code{numberProbes(object)}:

      Integer vector indicating the number of probes (markers) for each
      range in \code{object}.  Equivalent to \code{coverage2}.

    }

    \item{}{
      \code{state(object)}:

      Accessor for the \code{elementMetadata} column 'state', when
      applicable. State is used to contain the index of the inferred
      copy number state for various \code{hmm} methods defined in the
      \pkg{VanillaICE}.

    }

  }
  }

  \seealso{\code{\linkS4class{GRanges}}}

  \examples{
	  library(IRanges)
	  library(GenomicRanges)
	  gr1 <- GRanges(seqnames = "chr2", ranges = IRanges(3, 6),
			 state=3L, numberProbes=100L)
	  ## convenience functions
	  state(gr1)
	  numberProbes(gr1)

	  gr2 <- GRanges(seqnames = c("chr1", "chr1"),
			 ranges = IRanges(c(7,13), width = 3),
			 state=c(2L, 2L), numberProbes=c(200L, 250L))
	  gr3 <- GRanges(seqnames = c("chr1", "chr2"),
			 ranges = IRanges(c(1, 4), c(3, 9)),
			 state=c(1L, 4L), numberProbes=c(300L, 350L))
	  ## Ranges organized by sample
	  grl <- GRangesList("sample1" = gr1, "sample2" = gr2, "sample3" = gr3)
	  sampleNames(grl) ## same as names(grl)
	  numberProbes(grl)
	  chromosome(grl)
	  state(grl)
	  gr <- stack(grl)
	  sampleNames(gr)
	  chromosome(gr)
	  state(gr)
  }

\keyword{methods}
