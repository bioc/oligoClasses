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
	validAnn <- supportedAnnotation()
	x %in% validAnn
}

supportedAnnotation <- function(){
	c("pd.mapping50k.hind240", "pd.mapping50k.xba240",
	  "pd.mapping50k.hind240,pd.mapping50k.xba240",	  
	  "pd.mapping250k.nsp",
	  "pd.mapping250k.sty",
	  "pd.mapping250k.nsp,pd.mapping250k.sty",
	  "pd.genomewidesnp.5",
	  "pd.genomewidesnp.6",
	  "genomewidesnp6",
	  "genomewidesnp5",
	  "human370v1c",
	  "human370quadv3c",
	  "human550v3b",
	  "human650v3a",
	  "human610quadv1b",
	  "human660quadv1a",
	  "human1mduov3b")	  
}
