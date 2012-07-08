setMethod("initialize", signature(.Object="GenomeAnnotatedDataFrame"),
	  function(.Object,
		   position=integer(),
		   isSnp=vector("logical", length(position)),
		   chromosome=vector("integer", length(position)),
		   row.names=NULL,
		   data,
		   varMetadata,
		   ...){
		  if(missing(data))
			  data <- data.frame(isSnp=isSnp, position=position, chromosome=chromosome, ..., row.names=row.names)
		  if(missing(varMetadata)){
			  nms <- names(list(...))
			  varMetadata <- data.frame(labelDescription=c("SNP indicator", "physical position", "chromosome", nms))
		  }
		  .Object <- callNextMethod(.Object, data=data, varMetadata=varMetadata)
	  })

setMethod("updateObject", signature(object="GenomeAnnotatedDataFrame"),
	  function(object, ..., verbose=FALSE){
		  ##as(object, "GenomeAnnotatedDataFrame")
		  ADF2GDF(object)
	 })

setMethod("coerce", signature(from="AnnotatedDataFrame", to="GenomeAnnotatedDataFrame"),
	  function(from, to){
		  new("GenomeAnnotatedDataFrame",
		      isSnp=as.logical(from$isSnp),
		      position=as.integer(from$position),
		      chromosome=as.integer(from$chromosome),
		      row.names=featureNames(from))
	  })

ADF2GDF <- function(object){
	new("GenomeAnnotatedDataFrame",
	    isSnp=as.logical(object$isSnp),
	    position=as.integer(object$position),
	    chromosome=as.integer(object$chromosome),
	    row.names=featureNames(object))
}


isValidGenomeAnnotatedDataFrame <- function(object){
		    if(!all(c("isSnp", "position", "chromosome") %in% varLabels(object)))
			    return("'isSnp', 'position', and 'chromosome' are required varLabels of the AnnotatedDataFrame for features")
		    if(!is(chromosome(object), "integer") || !is(position(object), "integer"))
			    return("chromosome and position must be integers. See function 'chromosome2integer'")
		    if(!is.logical(isSnp(object))){
			    return("isSnp must be logical")
		    }
		    TRUE
}

setValidity("GenomeAnnotatedDataFrame",
	    function(object){
		    msg <- isValidGenomeAnnotatedDataFrame(object)
		    if(is.character(msg)) return(msg)
	    })

GenomeAnnotatedDataFrameFromMatrix <- function(object,
					       annotationPkg,
					       genome,
					       ...){
	dims <- dim(object)
	if (is.null(dims) || all(dims==0)){
		object <- GenomeAnnotatedDataFrameFrom(NULL, ...)
	} else {
		nms <- rownames(object)
		if(!is.null(nms)){
			if(any(is.na(nms))) warning("NA's in rownames")
		}
		if(length(annotationPkg)==0 | is.null(nms)){
			n <- dims[1]
			data <- data.frame(position=integer(n),
					   chromosome=integer(n),
					   isSnp=logical(n),
					   row.names=nms)
			object <- new("GenomeAnnotatedDataFrame", data=data)
		} else {
			stopifnot(isSupportedAnnotation(annotationPkg))
			is.pd <- isPdAnnotationPkg(annotationPkg)
			if(is.pd){
				object <- addFeatureAnnotation.pd2(annotationPkg, featureNames=rownames(object), genome=genome,...)
			} else {
				object <- addFeatureAnnotation.crlmm2(annotationPkg, featureNames=rownames(object), genome=genome, ...)
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
        GenomeAnnotatedDataFrameFrom(object[[eltNames[1]]], annotationPkg, ...)
}

setMethod("GenomeAnnotatedDataFrameFrom",
	  signature(object="ff_or_matrix"),
	  function(object, annotationPkg, genome="hg19", ...){
		  GenomeAnnotatedDataFrameFromMatrix(object=object, annotationPkg=annotationPkg, genome=genome, ...)
})

setMethod("GenomeAnnotatedDataFrameFrom",
	  signature(object="NULL"),
	  function(object, annotationPkg, genome, ...){
		  GenomeAnnotatedDataFrameFromNULL(object)
})

setMethod("GenomeAnnotatedDataFrameFrom",
	  signature(object="AssayData"),
	  function(object, annotationPkg, genome="hg19", ...){
		  GenomeAnnotatedDataFrameFromAssayData(object=object, annotationPkg=annotationPkg, genome=genome, ...)
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
				    value="logical"),
		 function(object, value){
			 Biobase::fData(object)$isSnp <-  value
			 object
		 })
setReplaceMethod("position", signature(object="GenomeAnnotatedDataFrame",
				       value="integer"),
		 function(object, value){
			 Biobase::fData(object)$position <-  value
			 object
		 })

isPdAnnotationPkg <- function(object) length(grep("pd.", object)) >= 1

addFeatureAnnotation.pd2 <- function(annotation, featureNames, genome){
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
		gb <- tolower(genomeBuild(annoObject))
		if(gb != genome){
			stop(paste("Genome build in package ",
				   annotation, " is ", gb, ", but build ",
				   genome, " is requested.", sep=""))
		}
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
			   isSnp=as.logical(fD[, "isSnp"]),
			   chromosome=as.integer(fD[, "chromosome"]),
			   position=as.integer(fD[, "position"]),
			   row.names=featureNames)
	return(featureData)
}




addFeatureAnnotation.crlmm2 <- function(object, featureNames, genome="", ...){
	nm <- grep("Crlmm", object, ignore.case=TRUE)
	if(length(nm) == 0){
		pkgname <- paste(object, "Crlmm", sep="")
	} else pkgname <- object
	path <- system.file("extdata", package=pkgname)
	if(path=="") stop("Are you sure ", pkgname, " is installed?")

	## Most of our annotation packages have only hg19 build
	snpBuilds <- list.files(path, pattern="snpProbes")
	multiple.builds <-  length(grep("_hg1[89].rda", snpBuilds)) > 1
	if(!multiple.builds & length(snpBuilds) > 1) snpBuilds <- "snpProbes.rda"
	hgbuild <- length(grep("_hg1[89].rda", snpBuilds)) > 0
	if(multiple.builds){
		loader(paste("cnProbes_", genome, ".rda", sep=""), pkgname=pkgname, envir=.oligoClassesPkgEnv)
		cnProbes <- get("cnProbes", envir=.oligoClassesPkgEnv)
		loader(paste("snpProbes_", genome, ".rda", sep=""), pkgname=pkgname, envir=.oligoClassesPkgEnv)
		snpProbes <- get("snpProbes", envir=.oligoClassesPkgEnv)
	} else {
		if(genome=="hg19" & !hgbuild) genome <- ""
		if(genome=="hg18" & !hgbuild) stop("hg18 build requested but only hg19 build is available")
		if(genome %in% c("hg18","hg19")) genome <- paste("_", genome, sep="")
		requestedBuild <- paste("snpProbes", genome, ".rda", sep="")
		match.arg(requestedBuild, snpBuilds)
		loader(requestedBuild, pkgname=pkgname, envir=.oligoClassesPkgEnv)

		cnBuilds <- list.files(path, pattern="cnProbes")
		cnRequestedBuild <- paste("cnProbes", genome, ".rda", sep="")
		match.arg(cnRequestedBuild, cnBuilds)
		loader(cnRequestedBuild, pkgname=pkgname, envir=.oligoClassesPkgEnv)
		cnProbes <- get("cnProbes", envir=.oligoClassesPkgEnv)
		snpProbes <- get("snpProbes", envir=.oligoClassesPkgEnv)
	}
	snpProbes <- snpProbes[rownames(snpProbes) %in% featureNames, , drop=FALSE]
	cnProbes <- cnProbes[rownames(cnProbes) %in% featureNames, , drop=FALSE]
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
	    isSnp=as.logical(isSnp[I]),
	    position=as.integer(position[I]),
	    chromosome=as.integer(chrom[I]),
	    row.names=featureNames)
}

cleancdfname <- function(x) strsplit(x, "Crlmm")[[1]][[1]]

isSupportedAnnotation <- function(x){
	validAnn <- annotationPackages()
	validAnn <- validAnn[-grep(",", validAnn)]
	stripCrlmm <- sapply(validAnn, cleancdfname)
	validAnn <- unique(c(validAnn, stripCrlmm))
	x <- strsplit(x, ",")[[1]]
	for(i in seq_along(x)) match.arg(x[i], validAnn)
	return(TRUE)
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
	  "humanomni1quadv1bCrlmm",
	  "gw6crlmm")
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
		istrue <- isSupportedAnnotation(x)
	} else {
		istrue <- sapply(x, isSupportedAnnotation)
		istrue <- all(istrue)
	}
	return(istrue)
}

setReplaceMethod("position", signature(object="oligoSnpSet", value="integer"),
		 function(object, value){
			 position(featureData(object)) <- value
			 object
		 })

setReplaceMethod("position", signature(object="GenomeAnnotatedDataFrame", value="integer"),
		 function(object, value){
			 object$position <- value
			 object
		 })
