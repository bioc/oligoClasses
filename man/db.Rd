\name{db}
\alias{db}
\alias{db-methods}
\alias{db,FeatureSet-method}
\alias{db,DBPDInfo-method}
\alias{db,SnpSet-method}
\alias{db,SnpCnvQSet-method}
\alias{db,SnpQSet-method}
\title{Get the connection to the SQLite Database}
\description{
  This function will return the SQLite connection to the database
  associated to objects used in oligo.
}
\usage{
db(object)
}
\arguments{
  \item{object}{Object of valid class. See methods.}
}
\section{Methods}{
\describe{
\item{object = "FeatureSet"}{object of class FeatureSet}
\item{object = "SnpCallSet"}{object of class SnpCallSet}
\item{object = "DBPDInfo"}{object of class DBPDInfo}
\item{object = "SnpLevelSet"}{object of class SnpLevelSet}
}}
\value{
  SQLite connection.
}
\author{Benilton Carvalho}
\examples{
## db(object)
}
\keyword{methods}
