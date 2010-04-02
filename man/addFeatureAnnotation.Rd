\name{addFeatureAnnotation}
\alias{addFeatureAnnotation}
\title{Add genomic annotation (chromosome, position) for several SNP platforms.}
\description{Adds chromosome, position, and an indicator for whether the locus is polymorphic.}
\usage{addFeatureAnnotation(object)}
\arguments{  \item{object}{ An object extending the eSet class.}}
\value{ An \code{AnnotatedDataFrame}.}
\author{R. Scharpf}
\examples{
	if(require(pd.genomewidesnp.6)){
		conn <- db(pd.genomewidesnp.6)
		dbListTables(conn)
		dbListFields(conn, "featureSet")
		## get 5 snp identifiers
		##sql <- "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' LIMIT 5"
		sql <- "SELECT man_fsetid FROM featureSet LIMIT 5"
		ids <- dbGetQuery(conn, sql)[[1]]
		A <- B <- matrix(rnorm(25), 5, 5, dimnames=list(ids, LETTERS[1:5])) 
		obj <- new("AlleleSet",
			   alleleA=A,
			   alleleB=B,
			   annotation="pd.genomewidesnp.6")
		featureData(obj) <- addFeatureAnnotation(obj)
		fData(obj)

		##check against annotation package
		##sql <- "SELECT man_fsetid, chrom, physical_pos FROM featureSet WHERE man_fsetid LIKE 'SNP%' LIMIT 5"
		##dbGetQuery(conn, sql)
	}
	if(require(genomewidesnp6Crlmm)){
		##alternatively, could use the Crlmm annotation package
		obj2 <- new("AlleleSet",
			   alleleA=A,
			   alleleB=B,
			   annotation="genomewidesnp6")
		featureData(obj2) <- addFeatureAnnotation(obj2)
		fData(obj2)
	}
}
\keyword{manip}

