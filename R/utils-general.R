checkAnnotation <- function(x){
	if(length(x) == 0) return(FALSE)
	x <- strsplit(x, ",")[[1]]
	if(length(x) == 1){
		return(isSupportedAnnotation(x))
	}
	if(length(x) == 2)
		return(all(sapply(x, isSupportedAnnotation)))
}

loader <- function(theFile, envir, pkgname){
	theFile <- file.path(system.file(package=pkgname),
			     "extdata", theFile)
	if (!file.exists(theFile))
		stop("File ", theFile, " does not exist in ", pkgname)
	load(theFile, envir=envir)
}

requireAnnotation <- function(pkgname, lib=.libPaths()[1], verbose=TRUE){
  stopifnot(is.character(pkgname), !missing(pkgname))
  status <- require(pkgname, character.only=TRUE, quietly=!verbose)
  if (!status)
    status <- pdPkgFromBioC(pkgname, lib=lib, verbose=verbose)
  status
}

isSupportedAnnotation <- function(x){
	validAnn <- annotationPackages()
	x %in% validAnn
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



## Package Downloader/Installer
## returns TRUE if installed and FALSE otherwise
pdPkgFromBioC <- function(pkgname, lib=.libPaths()[1], verbose=TRUE) {
  if (length(lib) > 1) {
    warning("Ignoring all but first element of argument lib.")
    lib <- lib[1]
  }

  if (verbose){
    message("Attempting to obtain '", pkgname, "' from BioConductor website.")
    message("Checking to see if your internet connection works...")
  }

  if (testBioCConnection()) {
    ## Check for file permissions
    if (file.access(lib, mode=0) < 0){
      if (verbose) message("Directory '", lib, "' does not seem to exist.")
      return(FALSE)
    }

    if (file.access(lib, mode=2) < 0){
      if (verbose) message("You do not have write access to '", lib, "'.")
      return(FALSE)
    }

    biocContribUrl <- sapply(biocReposList(), contrib.url)
    biocPkgs <- available.packages(biocContribUrl)
    if (! pkgname %in% biocPkgs[, "Package"]) {
      if (verbose)
        message("Package '", pkgname, "' was not found in the BioConductor repository.\n",
                "The 'pdInfoBuilder' package can often be used in situations like this.")
      return(FALSE)
    } else {
      install.packages(pkgname, lib=lib,
                       repos=biocReposList(),
                       dependencies=TRUE)
      status <- require(pkgname, character.only=TRUE, quietly=!verbose)
      if (status){
        return(TRUE)
      }else{
        if (verbose)
          message("There was a problem during download or installation.\n",
                  "Package '", pkgname, "' cannot be loaded. Please, try again.")
        return(FALSE)
      }
    }
  } else {
    if (verbose)
      message("Could not access the Bioconductor repository.\n",
              "Please check your internet connection.")
    return(FALSE)
  }
}



list.celfiles <-   function(..., listGzipped=FALSE){
    files <- list.files(...)
    if (listGzipped){
      return(files[grep("\\.[cC][eE][lL]\\.[gG][zZ]$|\\.[cC][eE][lL]$", files)])
    }else{
      return(files[grep("\\.[cC][eE][lL]$", files)])
    }
}

celfileDate <- function(filename) {
	h <- affyio::read.celfile.header(filename, info="full")
	date <- grep("/", strsplit(h$DatHeader, " ")[[1]], value=TRUE)
	if(length(date) < 1){
		##try something else
		results <- h$ScanDate
	} else{
		date <- strsplit(date, split="/")[[1]]
		CC <- ifelse(substr(date[3],1,1)=="9", "19", "20")
		results <- as.character(as.Date(paste(paste(CC, date[3], sep=""), date[1],
						      date[2], sep="-")))
	}
	results
}

addFeatureAnnotation <- function(object){
	if(length(grep("pd.", annotation(object))) >= 1){
		fD <- addFeatureAnnotation.pd(object)
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
	featureData <- new("AnnotatedDataFrame", data=fD,
			   varMetadata=data.frame(labelDescription=colnames(fD)))
	##Figure out how to add an indicator for SNP/CN probe
	return(featureData)
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

featureDataFrom <- function(annotationPackage){
	stopifnot(isSupportedAnnotation(annotationPackage))
	pkgname <- paste(cdfName, "Crlmm", sep="")
	path <- system.file("extdata", package=pkgname)
	loader("cnProbes.rda", pkgname=pkgname, envir=.oligoClassesPkgEnv)
	cnProbes <- get("cnProbes", envir=.oligoClassesPkgEnv)
	loader("snpProbes.rda", pkgname=pkgname, envir=.oligoClassesPkgEnv)
	snpProbes <- get("snpProbes", envir=.oligoClassesPkgEnv)
	fns <- c(rownames(snpProbes), rownames(cnProbes))
	isSnp <- c(rep(1L, nrow(snpProbes)), rep(0L, nrow(cnProbes)))
	positions <- as.integer(c(snpProbes[, "position"], cnProbes[, "position"]))
	chroms <- c(snpProbes[, "chrom"], cnProbes[, "chrom"])
	chroms <- chromosome2integer(chroms)
	tmp.fd <- cbind(chroms,  positions, isSnp)
	rownames(tmp.fd) <- fns
	tmp.fd <- tmp.fd[order(tmp.fd[, "chroms"], tmp.fd[, "positions"]), ]
	colnames(tmp.fd) <- c("chromosome", "position", "isSnp")
	featureData <- new("AnnotatedDataFrame", data=data.frame(tmp.fd), varMetadata=data.frame(labelDescription=colnames(tmp.fd)))
	return(featureData)
}

addFeatureAnnotation.crlmm <- function(object, ...){
	##if(missing(CHR)) stop("Must specificy chromosome")
	##message("Adding required feature annotation (chromosome, position, isSnp) to featureData slot")
	cdfName <- annotation(object)
	pkgname <- paste(cdfName, "Crlmm", sep="")
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
	fD <- new("AnnotatedDataFrame", data=fd, varMetadata=data.frame(labelDescription=colnames(fd), row.names=colnames(fd)))
	return(fD)
}

## a bar that I like to use when sending messages to the user
getBar <- function(width=getOption("width"))
  paste(rep("=", width), collapse="")

isPackageLoaded <- function(pkg){
  stopifnot(is.character(pkg))
  pkg <- paste("package:", pkg, sep="")
  pkg %in% search()
}


checkExists <- function(.name, .path=".", .FUN, .FUN2, .save.it=TRUE, .load.it, ...){
	##default of load.it depends on whether the object exists in .GlobalEnv
	if(exists(.name)){
		message("Exists in .GlobalEnv")
		if(missing(.load.it)){
			message(".load.it is missing. Setting .load.it to FALSE")
			.load.it <- FALSE
		}
		if(.load.it){
			fname <- file.path(.path, paste(.name, ".rda", sep=""))
			if(file.exists(fname)){
				message(".load.it is TRUE")
				message("Loading ", fname)
				load(fname)
				if(!exists(".object")) .object <- get(.name)
				return(.object)
			} else {
				message(fname, " does not exist")
				message("Running .FUN")
				.object <- .FUN(...)
				if(.save.it) {
					message("Saving ", fname)
					save(.object, file=fname)
				}
				return(.object)
			}
		} else {
			message(".load.it is FALSE. Nothing to do")
			.object <- get(.name)
			return(.object)
		}
	} else{
		message(.name, " does not exist in .GlobalEnv")
		fname <- file.path(.path, paste(.name, ".rda", sep=""))
		if(file.exists(fname)){
			message(fname, " exists")
			if(missing(.load.it)){
				message(".load.it is missing. Setting .load.it to TRUE")
				.load.it <- TRUE
			}
			if(.load.it){
				message("Loading ", fname)
				.tmp <- ls()
				load(fname)
				if(!exists(".object")) .object <- tryCatch(get(.name), error=function(e) NULL)
				##extremely ad-hoc
				if(is.null(.object)) .object <- get(ls()[!(ls() %in% .tmp) & !(ls() %in% c(".object", ".tmp"))])
				return(.object)
			} else {
				message(".load.it is FALSE.  Running .FUN")
				.object <- .FUN(...)
				if(.save.it) {
					message("Saving ", fname)
					save(.object, file=fname)
				}
				return(.object)
			}
		} else {
			message(fname, " does not exist. Running .FUN")
			.object <- .FUN(...)
			if(.save.it) {
				message("Saving ", fname)
				save(.object, file=fname)
			}
			return(.object)
		}
	}
}



