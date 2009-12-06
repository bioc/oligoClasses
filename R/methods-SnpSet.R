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

setMethod("callsConfidence", "SnpSet", function(object) confs(object))
setReplaceMethod("callsConfidence", signature(object="SnpSet", value="matrix"),
                 function(object, value) confs(object) <- value)

setMethod("isSnp", "SnpSet", function(object) {
	labels <- fvarLabels(object)
	if("isSnp" %in% labels){
		res <- fData(object)[, "isSnp"]
	} else{
		res <- as.integer(featureNames(object) %in% snpNames(object))
	}
	return(res==1)
})


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
	require(annotation(object), character.only=TRUE)
	message("Adding required feature annotation (chromosome, position, isSnp) to featureData slot")
	fs <- featureNames(object)
	tmp <- paste("('", paste(fs, collapse="', '"), "')", sep="")
	fD <- matrix(integer(), length(fs), 3)
	rownames(fD) <- fs
	colnames(fD) <- c("chromosome", "position", "isSnp")
	sql <- paste("SELECT man_fsetid, chrom, physical_pos FROM featureSet WHERE man_fsetid IN ", tmp)
	##Check if two objects have been combined
	pkgs <- strsplit(annotation(object), ",")[[1]]
	snps <- snp.index <- nps <- np.index <- vector("list", length(pkgs))
	for(i in seq(along=pkgs)){
		annotation(object) <- pkgs[i]
                snps[[i]] <- dbGetQuery(db(object), sql)
		snp.index[[i]] <- match(snps[[i]]$man_fsetid, rownames(fD))

		if("featureSetCNV" %in% dbListTables(db(object))){
			sql <- paste("SELECT man_fsetid, chrom, chrom_start FROM featureSetCNV WHERE man_fsetid IN ", tmp)
			nps[[i]] <- dbGetQuery(db(object), sql)
			np.index[[i]] <- match(nps[[i]]$man_fsetid, rownames(fD))
		}
	}
	if(length(snps) > 1){
		snps <- do.call(rbind, snps)
		snp.index <- unlist(snp.index)
		if("featureSetCNV" %in% dbListTables(db(object))){		
			nps <- do.call(rbind, nps)
			np.index <- unlist(np.index)
		}
	} else {
		snps <- snps[[1]]
		snp.index <- snp.index[[1]]		
		if("featureSetCNV" %in% dbListTables(db(object))){		
			nps <- nps[[1]]
			np.index <- np.index[[1]]
		}
	}
	fD[snp.index, "isSnp"] <- as.integer(1)
	fD[snp.index, "chromosome"] <- chromosome2integer(snps$chrom)
	fD[snp.index, "position"] <- as.integer(snps$physical_pos)
	if("featureSetCNV" %in% dbListTables(db(object))){			
		fD[np.index, "isSnp"] <- as.integer(0)
		fD[np.index, "chromosome"] <- chromosome2integer(nps$chrom)
		fD[np.index, "position"] <- as.integer(nps$chrom_start)
	}
	fD <- cbind(fD, fData(object))
	featureData <- new("AnnotatedDataFrame", data=fD,
			   varMetadata=data.frame(labelDescription=colnames(fD)))
	##Figure out how to add an indicator for SNP/CN probe
	return(featureData)
}

loader <- function(theFile, envir, pkgname){
	theFile <- file.path(system.file(package=pkgname),
			     "extdata", theFile)
	if (!file.exists(theFile))
		stop("File ", theFile, " does not exist in ", pkgname)
	load(theFile, envir=envir)
}

chromosome2integer <- function(chrom){
	chrom[chrom == "X"] <- 23; chrom[chrom == "Y"] <- 24; chrom[chrom == "XY"] <- 25; chrom[chrom=="M" | chrom == "MT"] <- 26
	as.integer(chrom)
}

snpNames <- function(object){
	path <- system.file("extdata", package=paste(annotation(object), "Crlmm", sep=""))
	load(file.path(path, "snpProbes.rda"))
	snpProbes <- get("snpProbes")
	snps <- rownames(snpProbes)
	snps <- snps[snps %in% featureNames(object)]
	index <- match(snps, featureNames(object), nomatch=0)
	index <- index[index != 0]
	featureNames(object)[index]
}

addFeatureAnnotation.crlmm <- function(object, ...){
	##if(missing(CHR)) stop("Must specificy chromosome")
	message("Adding required feature annotation (chromosome, position, isSnp) to featureData slot")
	cdfName <- annotation(object)
	pkgname <- paste(cdfName, "Crlmm", sep="")	
	path <- system.file("extdata", package=pkgname)
	loader("cnProbes.rda", pkgname=pkgname, envir=.oligoClassesPkgEnv)
	cnProbes <- get("cnProbes", envir=.oligoClassesPkgEnv)
	loader("snpProbes.rda", pkgname=pkgname, envir=.oligoClassesPkgEnv)
	snpProbes <- get("snpProbes", envir=.oligoClassesPkgEnv)	
	##Feature Data
	isSnp <- rep(as.integer(0), nrow(object))
	snpIndex <- function(object){
		index <- match(snpNames(object), featureNames(object), nomatch=0)
		index[index != 0]
	}
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
