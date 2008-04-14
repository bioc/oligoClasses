\name{db}
\alias{db}
\alias{db-methods}
\alias{db,SnpCallSet-method}
\alias{db,DBPDInfo-method}
\title{Get the connection to the SQLite Database}
\description{
  This function will return the SQLite connection to the database
  associated to objects used in oligo.
}
\usage{
db(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of valid class. See methods.}
}
\section{Methods}{
\describe{
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
