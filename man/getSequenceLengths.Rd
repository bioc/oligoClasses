\name{getSequenceLengths}
\alias{getSequenceLengths}

\title{Load chromosome sequence lengths for UCSC genome build hg18 or hg19}
\description{Load chromosome sequence lengths for UCSC genome build hg18 or hg19}
\usage{
getSequenceLengths(build)
}
\arguments{
  \item{build}{
    character string: "hg18" or "hg19"
  }
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

