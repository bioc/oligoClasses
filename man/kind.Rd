\name{kind}
\alias{kind}
\alias{kind,AffySNPPDInfo-method}
\alias{kind,AffySNPCNVPDInfo-method}
\alias{kind,ExpressionPDInfo-method}
\alias{kind,TilingPDInfo-method}
\alias{kind,AffyExpressionPDInfo-method}
\alias{kind,AffyGenePDInfo-method}
\alias{kind,FeatureSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Array type}
\description{
  Retrieves the array type.
}
\usage{
kind(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{FeatureSet} object}
}
\value{
  String: "Expression", "Exon", "SNP" or "Tiling"
}
\examples{
if (require(pd.mapping50k.xba240)){
 data(sfsExample)
 annotation(sfsExample) <- "pd.mapping50k.xba240"
 kind(sfsExample)
}
}

\keyword{manip}
