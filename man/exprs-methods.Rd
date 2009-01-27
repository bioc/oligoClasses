\name{exprs-methods}
\docType{methods}
\alias{exprs-methods}
%\alias{exprs,ExpressionSet-method}
\alias{exprs,FeatureSet-method}
\alias{se.exprs,FeatureSet-method}
%\alias{exprs,SnpSet-method}
\title{Accessor for the 'exprs' slot}
\description{
  Accessor for the 'exprs'/'se.exprs' slot of FeatureSet-like objects
}
\section{Methods}{
\describe{
\item{object = "ExpressionSet"}{Expression matrix for objects of this
  class. Usually results of preprocessing algorithms, like RMA.}
\item{object = "FeatureSet"}{General container 'exprs' inherited from eSet}
\item{object = "SnpSet"}{General container 'exprs' inherited from eSet,
  not yet used.}
}}
\keyword{methods}
