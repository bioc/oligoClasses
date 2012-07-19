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

    biocContribUrl <- sapply(biocinstallRepos(), contrib.url)
    biocPkgs <- available.packages(biocContribUrl)
    if (! pkgname %in% biocPkgs[, "Package"]) {
      if (verbose)
        message("Package '", pkgname, "' was not found in the BioConductor repository.\n",
                "The 'pdInfoBuilder' package can often be used in situations like this.")
      return(FALSE)
    } else {
      install.packages(pkgname, lib=lib,
                       repos=biocinstallRepos(),
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

celfileName <- function(object){
	if(!is(object, "CNSet")) stop("object must be CNSet")
	vl <- varLabels(protocolData(object))
	dirnames <- rep(object@datadir[[1]], object@datadir[[2]])
	if("filename" %in% vl){
		fns <- file.path(dirnames, protocolData(object)$filename)
	} else {
		fns <- file.path(dirnames, basename(sampleNames(object)))
	}
	if(!all(file.exists(fns))) stop("not all files exist")
	fns
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
				message("Running ", .FUN)
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

library2 <- function(...){
	suppressPackageStartupMessages(library(...))
}



