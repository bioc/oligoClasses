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
}
\keyword{attribute}

