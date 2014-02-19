\name{kind}
\alias{kind}
\alias{kind,AffySNPPDInfo-method}
\alias{kind,AffySNPCNVPDInfo-method}
\alias{kind,AffyExonPDInfo-method}
\alias{kind,ExpressionPDInfo-method}
\alias{kind,TilingPDInfo-method}
\alias{kind,AffyExpressionPDInfo-method}
\alias{kind,AffyGenePDInfo-method}
\alias{kind,AffyHTAPDInfo-method}
\alias{kind,FeatureSet-method}

\title{Array type}
\description{
  Retrieves the array type.
}
\usage{
kind(object)
}

\arguments{
  \item{object}{\code{FeatureSet} or \code{DBPDInfo} object}
}
\value{
  String: "Expression", "Exon", "SNP" or "Tiling"
}
\examples{
if (require(pd.mapping50k.xba240)){
 data(sfsExample)
 Biobase::annotation(sfsExample) <- "pd.mapping50k.xba240"
 kind(sfsExample)
}
}

\keyword{manip}
