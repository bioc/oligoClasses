setMethod("initialize", signature(.Object="GenomeAnnotatedDataFrame"),
	  function(.Object,
		   position=integer(),
		   isSnp=vector("integer", length(position)),
		   chromosome=vector("integer", length(position)),
		   row.names=NULL,
		   data=data.frame(isSnp=isSnp, position=position, chromosome=chromosome, row.names=row.names),
		   varMetadata=data.frame(labelDescription=c("SNP indicator", "physical position", "chromosome")),
		   ...){
		  rownames(varMetadata) <- c("isSnp", "position", "chromosome")
		  .Object <- callNextMethod(.Object, data=data, varMetadata=varMetadata, ...)
	  })

GenomeAnnotatedDataFrameFromMatrix <- function(object,
					    annotationPkg,
					    ...){
	dims <- dim(object)
	if (is.null(dims) || all(dims==0)){
		object <- GenomeAnnotatedDataFrameFrom(NULL, ...)
	} else {
		nms <- rownames(object)
		if(length(annotationPkg)==0){
			n <- dims[1]
			data <- data.frame(position=integer(n),
					   chromosome=integer(n),
					   isSnp=integer(n),
					   row.names=nms)
			object <- new("GenomeAnnotatedDataFrame", data=data)
		} else {
			stopifnot(isSupportedAnnotation(annotationPkg))
			is.pd <- isPdAnnotationPkg(annotationPkg)
			if(is.pd){
				object <- addFeatureAnnotation.pd2(annotationPkg, rownames(object))
			} else {
				object <- addFeatureAnnotation.crlmm2(annotationPkg, rownames(object))
			}
		}
	}
	return(object)
}

GenomeAnnotatedDataFrameFromNULL <- function(object, byrow=TRUE, ...) {
	new("GenomeAnnotatedDataFrame")
}

GenomeAnnotatedDataFrameFromAssayData <- function(object, annotationPkg, ...) {
    eltNames <-
        if (is(object, "environment")) ls(object)
        else names(object)
    if (length(eltNames)==0)
	    GenomeAnnotatedDataFrameFrom(NULL)
    else
        GenomeAnnotatedDataFrameFrom(object[[eltNames[1]]], annotationPkg)
}

setMethod("GenomeAnnotatedDataFrameFrom",
	  signature(object="matrix"),
	  function(object, annotationPkg){
		  GenomeAnnotatedDataFrameFromMatrix(object, annotationPkg)
})

setMethod("GenomeAnnotatedDataFrameFrom",
	  signature(object="NULL"),
	  function(object){
		  GenomeAnnotatedDataFrameFromNULL(object)
})

setMethod("GenomeAnnotatedDataFrameFrom",
	  signature(object="AssayData"),
	  function(object, annotationPkg){
		  GenomeAnnotatedDataFrameFromAssayData(object, annotationPkg)
})

setMethod("isSnp", signature(object="GenomeAnnotatedDataFrame"),
	  function(object) object$isSnp)
setMethod("position", signature(object="GenomeAnnotatedDataFrame"),
	  function(object, ...) object$position)
setMethod("chromosome", signature(object="GenomeAnnotatedDataFrame"),
	  function(object, ...) object$chromosome)
setReplaceMethod("chromosome", signature(object="GenomeAnnotatedDataFrame",
					 value="integer"),
		 function(object, value){
			 fData(object)$chromosome <-  value
			 object
		 })
setReplaceMethod("isSnp", signature(object="GenomeAnnotatedDataFrame",
				    value="integer"),
		 function(object, value){
			 fData(object)$isSnp <-  value
			 object
		 })
setReplaceMethod("position", signature(object="GenomeAnnotatedDataFrame",
				       value="integer"),
		 function(object, value){
			 fData(object)$position <-  value
			 object
		 })

isPdAnnotationPkg <- function(object) length(grep("pd.", object)) >= 1

addFeatureAnnotation <- function(object){
	if(length(grep("pd.", annotation(object))) >= 1){
		fD <- addFeatureAnnotation.pd(object)
	} else {
		fD <- addFeatureAnnotation.crlmm(object)
	}
	return(fD)
}

addFeatureAnnotation2 <- function(object){
	if(length(grep("pd.", annotation(object))) >= 1){
	} else {
		fD <- addFeatureAnnotation.crlmm(object)
	}
	return(fD)
}

addFeatureAnnotation.pd <- function(object){
	##message("Adding required feature annotation (chromosome, position, isSnp) to featureData slot")
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
		require(annotation(object), character.only=TRUE)
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
	jj <- match(c("chromosome", "position", "isSnp"), fvarLabels(object))
	jj <- jj[!is.na(jj)]
	if(length(jj) > 0){
		fD <- cbind(fD, fData(object)[, -jj, drop=FALSE])
	} else fD <- cbind(fD, fData(object))
	featureData <- new("GenomeAnnotatedDataFrame",
			   isSnp=fD[, "isSnp"],
			   chromosome=fD[, "chromosome"],
			   position=fD[, "position"])
	##Figure out how to add an indicator for SNP/CN probe
	return(featureData)
}

addFeatureAnnotation.pd2 <- function(annotation, featureNames){
	##message("Adding required feature annotation (chromosome, position, isSnp) to featureData slot")
	fs <- featureNames
	tmp <- paste("('", paste(fs, collapse="', '"), "')", sep="")
	fD <- matrix(integer(), length(fs), 3)
	rownames(fD) <- fs
	colnames(fD) <- c("chromosome", "position", "isSnp")
	sql <- paste("SELECT man_fsetid, chrom, physical_pos FROM featureSet WHERE man_fsetid IN ", tmp)
	##Check if two objects have been combined
	pkgs <- strsplit(annotation, ",")[[1]]
	snps <- snp.index <- nps <- np.index <- vector("list", length(pkgs))
	for(i in seq(along=pkgs)){
		annotation <- pkgs[i]
		requireAnnotation(annotation)
		annoObject <- get(annotation)
                snps[[i]] <- dbGetQuery(annoObject@getdb(), sql)
		snp.index[[i]] <- match(snps[[i]]$man_fsetid, rownames(fD))
		if("featureSetCNV" %in% dbListTables(annoObject@getdb())){
			sql <- paste("SELECT man_fsetid, chrom, chrom_start FROM featureSetCNV WHERE man_fsetid IN ", tmp)
			nps[[i]] <- dbGetQuery(annoObject@getdb(), sql)
			np.index[[i]] <- match(nps[[i]]$man_fsetid, rownames(fD))
		}
	}
	if(length(snps) > 1){
		snps <- do.call(rbind, snps)
		snp.index <- unlist(snp.index)
		if("featureSetCNV" %in% dbListTables(annoObject@getdb())){
			nps <- do.call(rbind, nps)
			np.index <- unlist(np.index)
		}
	} else {
		snps <- snps[[1]]
		snp.index <- snp.index[[1]]
		if("featureSetCNV" %in% dbListTables(annoObject@getdb())){
			nps <- nps[[1]]
			np.index <- np.index[[1]]
		}
	}
	fD[snp.index, "isSnp"] <- as.integer(1)
	fD[snp.index, "chromosome"] <- chromosome2integer(snps$chrom)
	fD[snp.index, "position"] <- as.integer(snps$physical_pos)
	if("featureSetCNV" %in% dbListTables(annoObject@getdb())){
		fD[np.index, "isSnp"] <- as.integer(0)
		fD[np.index, "chromosome"] <- chromosome2integer(nps$chrom)
		fD[np.index, "position"] <- as.integer(nps$chrom_start)
	}
	featureData <- new("GenomeAnnotatedDataFrame",
			   isSnp=as.integer(fD[, "isSnp"]),
			   chromosome=as.integer(fD[, "chromosome"]),
			   position=as.integer(fD[, "position"]),
			   row.names=featureNames)
	return(featureData)
}

chromosome2integer <- function(chrom){
	chrom[chrom == "X"] <- 23; chrom[chrom == "Y"] <- 24; chrom[chrom == "XY"] <- 25; chrom[chrom=="M" | chrom == "MT" | chrom == "Mt"] <- 26
	as.integer(chrom)
}

featureDataFrom <- function(annotationPackage){
	cdfName <- strsplit(annotationPackage, "Crlmm")[[1]][[1]]
	pkgname <- paste(cdfName, "Crlmm", sep="")
	stopifnot(isSupportedAnnotation(pkgname))
	path <- system.file("extdata", package=pkgname)
	loader("cnProbes.rda", pkgname=pkgname, envir=.oligoClassesPkgEnv)
	cnProbes <- get("cnProbes", envir=.oligoClassesPkgEnv)
	loader("snpProbes.rda", pkgname=pkgname, envir=.oligoClassesPkgEnv)
	snpProbes <- get("snpProbes", envir=.oligoClassesPkgEnv)
	if("chr" %in% colnames(snpProbes)) {
		colnames(cnProbes) <- colnames(snpProbes) <- c("chrom", "position")
	}
	fns <- c(rownames(snpProbes), rownames(cnProbes))
	isSnp <- c(rep(1L, nrow(snpProbes)), rep(0L, nrow(cnProbes)))
	positions <- as.integer(c(snpProbes[, "position"], cnProbes[, "position"]))
	chroms <- c(snpProbes[, "chrom"], cnProbes[, "chrom"])
	chroms <- chromosome2integer(chroms)
	tmp.fd <- cbind(chroms,  positions, isSnp)
	rownames(tmp.fd) <- fns
	tmp.fd <- tmp.fd[order(tmp.fd[, "chroms"], tmp.fd[, "positions"]), ]
	colnames(tmp.fd) <- c("chromosome", "position", "isSnp")
	##featureData <- new("GenomeAnnotatedDataFrame", data=data.frame(tmp.fd), varMetadata=data.frame(labelDescription=colnames(tmp.fd)))
	featureData <- new("GenomeAnnotatedDataFrame",
			   isSnp=tmp.fd[, "isSnp"],
			   position=tmp.fd[, "position"],
			   chromosome=tmp.fd[, "chromosome"])
	return(featureData)
}

addFeatureAnnotation.crlmm <- function(object, ...){
	##if(missing(CHR)) stop("Must specificy chromosome")
	##message("Adding required feature annotation (chromosome, position, isSnp) to featureData slot")
	cdfName <- annotation(object)
	nm <- grep("Crlmm", cdfName)
	if(length(nm) == 0){
		pkgname <- paste(cdfName, "Crlmm", sep="")
	} else pkgname <- cdfName
	path <- system.file("extdata", package=pkgname)
	loader("cnProbes.rda", pkgname=pkgname, envir=.oligoClassesPkgEnv)
	cnProbes <- get("cnProbes", envir=.oligoClassesPkgEnv)
	loader("snpProbes.rda", pkgname=pkgname, envir=.oligoClassesPkgEnv)
	snpProbes <- get("snpProbes", envir=.oligoClassesPkgEnv)
	##Feature Data
	isSnp <- 1L-as.integer(featureNames(object) %in% rownames(cnProbes))
	names(isSnp) <- featureNames(object)
	if(any(isSnp)){
		snps <- featureNames(object)[isSnp == 1]
		position.snp <- snpProbes[match(snps, rownames(snpProbes)), "position"]
		names(position.snp) <- snps
		J <- grep("chr", colnames(snpProbes))
		chr.snp <- snpProbes[match(snps, rownames(snpProbes)), J]
	} else{
		warning("None of the featureNames in the object match SNP probes for the indicated annotation package.  Either the annotation package is misspecified, or the featureNames of the object are incorrect")
		message("The first 5 featureNames are ", featureNames(object)[1:5])
		message("The annotation for the object is ", annotation(object))
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
		warning("physical position not available for all featureNames")
		## Very dangerous with ff objects
		##object <- object[featureNames(object) %in% names(position), ]
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
		tmp.fd <- tmp.fd[, -grep("chromosome", colnames(tmp.fd)), drop=FALSE]
	if("position" %in% fvarLabels(object))
		tmp.fd <- tmp.fd[, -grep("position", colnames(tmp.fd)), drop=FALSE]
	if("isSnp" %in% fvarLabels(object))
		tmp.fd <- tmp.fd[, -grep("isSnp", colnames(tmp.fd)), drop=FALSE]
	rownames(tmp.fd) <- featureNames(object)
	tmp <- new("AnnotatedDataFrame",
		   data=tmp.fd,
		   varMetadata=data.frame(labelDescription=colnames(tmp.fd)))
	fd <- cbind(pData(tmp), fData(object))
	##fD <- new("AnnotatedDataFrame", data=fd, varMetadata=data.frame(labelDescription=colnames(fd), row.names=colnames(fd)))
	fD <- new("GenomeAnnotatedDataFrame",
		  isSnp=fd$isSnp,
		  position=fd$position,
		  chromosome=fd$chromosome)
	return(fD)
}

addFeatureAnnotation.crlmm2 <- function(object, featureNames, ...){
	##if(missing(CHR)) stop("Must specificy chromosome")
	##message("Adding required feature annotation (chromosome, position, isSnp) to featureData slot")
	##cdfName <- annotation(object)
	nm <- grep("Crlmm", object)
	if(length(nm) == 0){
		pkgname <- paste(object, "Crlmm", sep="")
	} else pkgname <- object
	path <- system.file("extdata", package=pkgname)
	loader("cnProbes.rda", pkgname=pkgname, envir=.oligoClassesPkgEnv)
	cnProbes <- get("cnProbes", envir=.oligoClassesPkgEnv)
	loader("snpProbes.rda", pkgname=pkgname, envir=.oligoClassesPkgEnv)
	snpProbes <- get("snpProbes", envir=.oligoClassesPkgEnv)

	snpProbes <- snpProbes[rownames(snpProbes) %in% featureNames, ]
	cnProbes <- cnProbes[rownames(snpProbes) %in% featureNames, ]
	##Feature Data
	isSnp <- 1L-as.integer(featureNames %in% rownames(cnProbes))
	names(isSnp) <- featureNames
	if(any(isSnp)){
		snps <- featureNames[isSnp == 1]
		index <- match(snps, rownames(snpProbes))
		position.snp <- snpProbes[index, "position"]
		names(position.snp) <- snps
		J <- grep("chr", colnames(snpProbes))
		chr.snp <- snpProbes[index, J]
	} else{
		warning("None of the featureNames in the object match SNP probes for the indicated annotation package.  Either the annotation package is misspecified, or the featureNames of the object are incorrect")
		##message("The first 5 featureNames are ", featureNames(object)[1:5])
		message("The annotation for the object is ", object)
		chr.snp <- position.snp <- integer()
	}
	if(any(!isSnp)){
		nps <- featureNames[isSnp == 0]
		index <- match(nps, rownames(cnProbes))
		position.np <- cnProbes[index, "position"]
		names(position.np) <- nps
		chr.np <- cnProbes[index, J]
	} else {
		chr.np <- position.np <- integer()
	}
	position <- c(position.snp, position.np)
	chrom <- c(chr.snp, chr.np)
	##We may not have annotation for all of the snps
	if(!all(featureNames %in% names(position))){
		warning("physical position not available for all featureNames")
	}
	ix <- match(featureNames, names(position))
	position <- position[ix]
	chrom <- chrom[ix]
	##require(SNPchip)
	chrom <- chromosome2integer(chrom)
	stopifnot(identical(names(position), featureNames))
	if(sum(duplicated(names(position))) > 0){
		warning("Removing rows with NA identifiers...")
		##RS: fix this
		I <- which(!is.na(names(position)))
	}  else I <- seq(along=names(position))
##	tmp.fd <- data.frame(cbind(chrom[I],
##				   position[I],
##				   isSnp[I]), row.names=featureNames)
##	colnames(tmp.fd) <- c("chromosome", "position", "isSnp")
	new("GenomeAnnotatedDataFrame",
	    isSnp=isSnp[I],
	    position=position[I],
	    chromosome=chrom[I],
	    row.names=featureNames)
}

isSupportedAnnotation <- function(x){
	validAnn <- annotationPackages()
	L <- grep(x, validAnn)>=1
	if(L==1) return(TRUE)
	if(L > 1){
		L <- grep(paste(x, "Crlmm", sep=""), validAnn)
		res <- if(L==1) TRUE else FALSE
		return(res)
	}
}

annotationPackages <- function(){
	c("pd.mapping50k.hind240", "pd.mapping50k.xba240",
	  "pd.mapping50k.hind240,pd.mapping50k.xba240",
	  "pd.mapping250k.nsp",
	  "pd.mapping250k.sty",
	  "pd.mapping250k.nsp,pd.mapping250k.sty",
	  "pd.genomewidesnp.5",
	  "pd.genomewidesnp.6",
	  "genomewidesnp6Crlmm",
	  "genomewidesnp5Crlmm",
	  "human370v1cCrlmm",
	  "human370quadv3cCrlmm",
	  "human550v3bCrlmm",
	  "human650v3aCrlmm",
	  "human610quadv1bCrlmm",
	  "human660quadv1aCrlmm",
	  "human1mduov3bCrlmm",
	  "humanomni1quadv1bCrlmm")
}

affyPlatforms <- function(){
  platforms <- c("pd.mapping50k.xba240",
                 "pd.mapping50k.hind240",
                 "pd.mapping250k.nsp",
                 "pd.mapping250k.sty",
                 "pd.genomewidesnp.5",
                 "pd.genomewidesnp.6")
  combined <- rep(NA, 2)
  combined[1] <- paste(sort(platforms[1:2]), collapse=",")
  combined[2] <- paste(sort(platforms[3:4]), collapse=",")
  platforms <- c(platforms, combined)
  platforms
}


checkAnnotation <- function(x){
	if(length(x) == 0) return(FALSE)
	x <- strsplit(x, ",")[[1]]
	if(length(x) == 1){
		return(isSupportedAnnotation(x))
	}
	if(length(x) == 2)
		return(all(sapply(x, isSupportedAnnotation)))
}