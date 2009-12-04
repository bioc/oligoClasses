setMethod("calls", "SnpSet", function(object) assayData(object)$call)
setReplaceMethod("calls", signature(object="SnpSet", value="matrix"), function(object, value) assayDataElementReplace(object, "call", value))
setMethod("confs", "SnpSet", function(object) {
	X <- assayData(object)$callProbability
	P <- 1-exp(-X/1000)
	return(P)
})
setReplaceMethod("confs", signature(object="SnpSet", value="matrix"),
		 function(object, value){
			 ##convert probability to integer
			 P <- value
			 dns <- dimnames(P)
			 X <- -1000*log(1-P)
			 X <- matrix(as.integer(X), nrow(X), ncol(X))
			 dimnames(X) <- dns
			 assayDataElementReplace(object, "callProbability", X)
		 })

setMethod("callsConfidence", "SnpSet", function(object)
          assayDataElement(object, "callProbability"))

setReplaceMethod("callsConfidence", signature(object="SnpSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "callProbability", value))


setMethod("db", "SnpSet",
          function(object) {
		  requireAnnotation(annotation(object)) || stop(paste(annotation(object), "package not available"))
		  get(annotation(object))@getdb()
	  })

addFeatureAnnotation <- function(object){
	if(length(grep("pd.", annotation(object))) >= 1){
		fD <- addFeatureAnnotation.pd(object)
	} else {
		fD <- addFeatureAnnotation.crlmm(object)
	}
	return(fD)
}

addFeatureAnnotation.pd <- function(object){
	require("RSQLite") || stop("RSQLite package not available")
	require(annotation(object), character.only=TRUE)
	message("Adding required feature annotation (chromosome, position, isSnp) to featureData slot")
	fs <- featureNames(object)
	tmp <- paste("'", paste(fs, collapse="', '"), "'", sep="")

	##extracts the snps
	##sql <- "SELECT man_fsetid, chrom, physical_pos FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
	sql <- paste("SELECT man_fsetid, chrom, physical_pos FROM featureSet WHERE man_fsetid IN ", tmp)
	conn <- db(annotation(object))	
	##Check if two objects have been combined
	pkgs <- strsplit(annotation(object), ",")[[1]]
	if(length(pkgs) > 1){
                object2 <- object1 <- object
                annotation(object1) <- pkgs[1]
                annotation(object2) <- pkgs[2]

                tmp1 <- dbGetQuery(db(object1), sql)
                tmp2 <- dbGetQuery(db(object2), sql)
                tmp <- rbind(tmp1, tmp2)
	} else {
                tmp <- dbGetQuery(db(object), sql)
	}
	idx <- match(fs, tmp[["man_fsetid"]])
	featureAnn <- tmp[idx, c("chrom", "physical_pos")]
	colnames(featureAnn) <- c("chromosome", "position")
	fD <- fData(object)
	fD2 <- cbind(fD, featureAnn)
	featureData <- new("AnnotatedDataFrame", data=fD2,
			   varMetadata=data.frame(labelDescription=colnames(featureData)))

	##Figure out how to add an indicator for SNP/CN probe
	return(featureData)
}

addFeatureAnnotation.crlmm <- function(object, ...){
	##if(missing(CHR)) stop("Must specificy chromosome")
	cdfName <- annotation(object)
	pkgname <- paste(cdfName, "Crlmm", sep="")	
	path <- system.file("extdata", package=pkgname)
	loader("cnProbes.rda", pkgname=pkgname, envir=.crlmmPkgEnv)
	cnProbes <- get("cnProbes", envir=.crlmmPkgEnv)
	loader("snpProbes.rda", pkgname=pkgname, envir=.crlmmPkgEnv)
	snpProbes <- get("snpProbes", envir=.crlmmPkgEnv)	

	##Feature Data
	isSnp <- rep(as.integer(0), nrow(object))
	isSnp[snpIndex(object)] <- as.integer(1)
	names(isSnp) <- featureNames(object)

	if(any(isSnp)){
		snps <- featureNames(object)[isSnp == 1]
		position.snp <- snpProbes[match(snps, rownames(snpProbes)), "position"]
		names(position.snp) <- snps

		J <- grep("chr", colnames(snpProbes))
		chr.snp <- snpProbes[match(snps, rownames(snpProbes)), J]		
	} else{
		chr.snp <- position.snp <- integer()
	}
	if(any(!isSnp)){
		nps <- featureNames(object)[isSnp == 0]
		position.np <- cnProbes[match(nps, rownames(cnProbes)), "position"]
		names(position.np) <- nps
		
		chr.np <- cnProbes[match(nps, rownames(cnProbes)), J]	
	} else {
		chr.np <- position.np <- integer()
	}
	position <- c(position.snp, position.np)
	chrom <- c(chr.snp, chr.np)

	##We may not have annotation for all of the snps
	if(!all(featureNames(object) %in% names(position))){
		message("Dropping loci for which physical position  is not available.")
		object <- object[featureNames(object) %in% names(position), ]
	}
	ix <- match(featureNames(object), names(position))
	position <- position[ix]
	chrom <- chrom[ix]
	##require(SNPchip)
	chrom <- chromosome2integer(chrom)

	stopifnot(identical(names(position), featureNames(object)))
	if(sum(duplicated(names(position))) > 0){
		warning("Removing rows with NA identifiers...")
		##RS: fix this
		I <- which(!is.na(names(position)))
	}  else I <- seq(along=names(position))
	tmp.fd <- data.frame(cbind(chrom[I],
				   position[I],
				   isSnp[I]))
	colnames(tmp.fd) <- c("chromosome", "position", "isSnp")
	if("chromosome" %in% fvarLabels(object))
		tmp.fd <- tmp.fd[, -1]
	if("position" %in% fvarLabels(object))
		tmp.fd <- tmp.fd[, -grep("position", colnames(tmp.fd)), drop=FALSE]
	if("isSnp" %in% fvarLabels(object))
		tmp.fd <- tmp.fd[, -grep("isSnp", colnames(tmp.fd)), drop=FALSE]
	rownames(tmp.fd) <- featureNames(object)
	tmp <- new("AnnotatedDataFrame",
		   data=tmp.fd,
		   varMetadata=data.frame(labelDescription=colnames(tmp.fd)))
	fd <- cbind(pData(tmp), fData(object))
	fD <- new("AnnotatedDataFrame", data=fd, varMetadata=data.frame(labelDescription=colnames(fd), row.names=colnames(fd)))
	return(fD)
}

setMethod("chromosome", "SnpSet",
	  function(object){
		  if(!("chromosome" %in% fvarLabels(object))){
			  stop("chromosome not in fvarLabels")
		  } 
		  return(featureData(object)$chromosome)
	  })

setReplaceMethod("chromosome", c("SnpSet", "ANY"),
		 function(object, value){
			 fData(object)$chromosome <-  value
			 object
		 })

setMethod("position", "SnpSet",
          function(object){
		  if(!("position" %in% fvarLabels(object))){
			  stop("position not in fvarLabels")
		  }
		  pos <- featureData(object)$position
          })

setMethod("combine", signature=signature(x="SnpSet", y="SnpSet"),
          function(x, y, ...){
		  ##Check that both x and y are valid objects
		  if(!validObject(x)) stop("x is not a valid object")
		  if(!validObject(y)) stop("y is not a valid object")
		  annot <- paste(sort(c(annotation(x), annotation(y))), collapse=",")
		  annotation(x) <- annotation(y) <- annot

		  if(class(x) != class(y)){
			  stop("objects must have the same class")
		  }            
		  if(storageMode(assayData(x)) != storageMode(assayData(y))){
			  stop("objects must have same storage mode for assayData")
		  }

		  fd <- combine(featureData(x), featureData(y))
		  pd <- combine(phenoData(x), phenoData(y))            
		  ad.x <- as.list(assayData(x))
		  ad.y <- as.list(assayData(y))
		  ad.xy <- mapply(rbind, ad.x, ad.y, SIMPLIFY=FALSE)
		  id.x <- match(rownames(ad.xy[[1]]), featureNames(fd))
		  ee <- combine(experimentData(x), experimentData(y))
		  assayData(x) <- ad.xy
		  storageMode(assayData(x)) <- storageMode(assayData(y))            
		  experimentData(x) <- ee
		  featureData(x) <- fd
		  phenoData(x) <- pd
		  x
          })
