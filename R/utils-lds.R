## utilities for Large Dataset Support
##
## Summary:
##   - is.ffmatrix: test specifically for ff_matrix
##   - ldStatus: TRUE if Large Dataset Support is available
##   - ldPath: where to save ff files
##   - createFF: creates an ff object setting path appropriately
##               (leaving out of fftempdir b/c parallel processes
##               access the object very easily)


initializeBigArray <- function(name=basename(tempfile()), dim=c(0L,0L,0L), vmode="integer", initdata=NA){
	if(isPackageLoaded("ff")){
		results <- createFF(name=name,
				    dim=dim,
				    vmode=vmode, initdata=initdata)
	} else {
		init <- switch(vmode,
			       integer=as.integer(initdata),
			       double=as.double(initdata),
			       character=as.character(initdata),
			       stop("Mode ", vmode, " not implemented for regular matrices"))
		results <- array(init, dim=dim)
	}
	return(results)
}


initializeBigMatrix <- function(name=basename(tempfile()), nr=0L, nc=0L, vmode="integer", initdata=NA){
  if(isPackageLoaded("ff") & (nr > 0 | nc > 0)){
    if(prod(nr, nc) > 2^31){
      ##Need multiple matrices
      ## -- use ffdf
      ## How many samples per ff object
      S <- floor(2^31/nr - 1)
      ## How many ff objects
      L <- ceiling(nc/S)
      name <- paste(name, 1:L, sep="_")
      resultsff <- vector("list", L)
      for(i in 1:(L-1)){  ## the Lth object may have fewer than nc columns
        resultsff[[i]] <- createFF(name=name[i],
                                   dim=c(nr, S),
                                   vmode=vmode, initdata=initdata)
      }
      ##the Lth element
      leftOver <- nc - ((L-1)*S)
      resultsff[[L]] <- createFF(name=name[L],
                                 dim=c(nr, leftOver),
                                 vmode=vmode, initdata=initdata)
      results <- do.call(ffdf, resultsff)
      rm(resultsff); gc()
    } else {
      results <- createFF(name=name,
                          dim=c(nr, nc),
                          vmode=vmode, initdata=initdata)
    }
   }  else {
    init <- switch(vmode,
                   integer=as.integer(initdata),
                   double=as.double(initdata),
                   character=as.character(initdata),
                   stop("Mode ", vmode, " not implemented for regular matrices"))
    results <- matrix(init, nr, nc)
}
  return(results)
}

initializeBigVector <- function(name=basename(tempfile()), n=0L, vmode="integer", initdata=NA){
  if(isPackageLoaded("ff")){
    results <- ff(initdata=initdata, vmode=vmode, length=n,
                  pattern=file.path(ldPath(), basename(name)))
  }  else {
    init <- switch(vmode,
                   integer=as.integer(initdata),
                    double=as.double(initdata),
                    character=as.character(initdata),
                    stop("Mode ", vmode, " not implemented for regular matrices"))
    results <- rep(init, n)
  }
  return(results)
}



createFF <- function(name, dim, vmode="double", initdata=NULL)
	ff(initdata=initdata, vmode=vmode, dim=dim, pattern=file.path(ldPath(), basename(name)))

## TODO: really safe to use :::?
setMethod("annotatedDataFrameFrom", "ff_matrix",
          Biobase:::annotatedDataFrameFromMatrix)

is.ffmatrix <- function(object)
  is(object, "ff_matrix")

isFF <- function(object){
	names <- ls(assayData(object))
	is(assayData(object)[[names[[1]]]], "ff") | is(assayData(object)[[names[[1]]]], "ffdf")
}

ldPath <- function(path){
  if (missing(path)){
    return(getOption("ldPath"))
  }else{
	  if(!is.character(path)) stop("path is not a character string")
##    stopifnot(is.character(path))
    options(ldPath=path)
  }
}

ldSetOptions <- function(nsamples=100, nprobesets=20000,
                         path=getwd(), verbose=FALSE){
  ocProbesets(nprobesets)
  ocSamples(nsamples)
  ldPath(path)
  ldStatus(verbose)
  TRUE
}

ldStatus <- function(verbose=FALSE){
  ld <- isPackageLoaded("ff")
  if (verbose){
    message(getBar())
    message("Large dataset support for 'oligo/crlmm': ", appendLF=FALSE)
    if (ld){
      message("Enabled")
      ns <- prettyNum(c(ocProbesets(), ocSamples()), big.mark=",")
      message("    - Probesets: ", ns[1])
      message("    - Samples..: ", ns[2])
      message("    - Path.....: ", ldPath())
    }else{
      message("Disabled")
      message("     - Load 'ff'")
    }
    message(getBar())
  }
  return(ld)
}

## does nothing if not an ff object
setMethod("open", "numeric", function(con, ...) return(NULL))
setMethod("open", "matrix", function(con, ...) return(NULL))
setMethod("open", "array", function(con, ...) return(NULL))
setMethod("close", "numeric", function(con, ...) return(NULL))
setMethod("close", "matrix", function(con, ...) return(NULL))
setMethod("close", "array", function(con, ...) return(NULL))
