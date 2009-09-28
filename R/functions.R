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
      
    biocContribUrl <- sapply(Biobase:::biocReposList(), contrib.url)
    biocPkgs <- available.packages(biocContribUrl)
    if (! pkgname %in% biocPkgs[, "Package"]) {
      if (verbose) 
        message("Package '", pkgname, "' was not found in the BioConductor repository.\n",
                "The 'pdInfoBuilder' package can often be used in situations like this.")
      return(FALSE)
    } else {
      install.packages(pkgname, lib=lib,
                       repos=Biobase:::biocReposList(),
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

requireAnnotation <- function(pkgname, lib=.libPaths()[1], verbose=TRUE){
  stopifnot(is.character(pkgname), !missing(pkgname))
  status <- require(pkgname, character.only=TRUE, quietly=!verbose)
  if (!status){
    status <- pdPkgFromBioC(pkgname, lib=lib, verbose=verbose)
  }
  status
}

